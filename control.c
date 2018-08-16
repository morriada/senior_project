/*
 * Initialization of System - Initialize peripherals and communication.
 * Copyright (C) 2018 Adam Morris <morriada@mail.gvsu.edu>
 * Copyright (C) 2018 Nicholas Borchardt <borcharn@mail.gvsu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// Include Project Libraries
#include "control.h"

// Declare variables
static volatile uint32_t freq[] = { 
	148000000,
	149000000,
	150000000,
	151000000
};

int i, r;
int auto_gain = 0;
int disable_dither = 0;
int i2c_repeater_on = 1;
int i2c_repeater_off = 0;
uint32_t i2c_addr = 0;
uint8_t bias_data = 0;
uint8_t bias_noise = 1;
uint8_t i2c_value = 0;
uint32_t sample_rate = 2400000;
uint32_t if_freq = 0;

void usage(void)
{
  printf(
  	"TODO:Usage Info\n"
	);
  exit(1);
}

void rtlsdr_setup(void)
{
	for(i = 0; i < NUM_SDRS; ++i)
	{
		sdrs[i].id = i;
		rtlsdr_dev_t *dev = sdrs[i].dev;
		r = rtlsdr_open(&dev, sdrs[i].id);
	  
  	// Set the sample rate of the rtl-sdr
  	rtlsdr_set_sample_rate(dev, sample_rate);
  	// Disable dithering
  	rtlsdr_set_dithering(dev, disable_dither);
  	// Set the IF frequency
  	rtlsdr_set_if_freq(dev, if_freq);
  	// Set the center frequency
  	rtlsdr_set_center_freq(dev, freq[i]);
  	// Set the tuner gain mode to automatic
  	rtlsdr_set_tuner_gain_mode(dev, auto_gain);
	}
}

void file_save(int sdr_num)
{
	// Save buffer to file for sdr_num
}

void * collect_t(void *ptr)
{
	// Cast thread pointer
	int * idx = (int *)ptr;
	// Declare thread variables
	int ret, blocksize, n_read;

	ret = n_read = 0;
	blocksize = sdrs[*idx].blocksize;
	ret = rtlsdr_read_sync(sdrs[*idx].dev, sdrs[*idx].buffer, blocksize, &n_read);
	if(ret < 0) {
		fprintf(stderr, "Runtime error: %d at %s:%d\n", ret, __FILE__, __LINE__);
	} else if(n_read < blocksize) {
		fprintf(stderr, "Short read %d: %d/%d\n", sdrs[*idx].id, n_read, blocksize);
	}

	file_save(*idx);

	return NULL;
}

void rtlsdr_calibration(void)
{
	for(i = 0; i < NUM_SDRS; ++i)
	{
		rtlsdr_dev_t *dev = sdrs[i].dev;

		// Set the bias tee by setting the gpio bit 0 to bias_off
  	rtlsdr_set_bias_tee(dev, bias_noise);
  	// Set rtlsdr repeater for the i2communication via RTL2838
  	rtlsdr_set_i2c_repeater(dev, i2c_repeater_on);
  	// Set register to the output
  	rtlsdr_i2c_write_reg(dev, i2c_addr, 0x03, 00);
  	// Set value to the register as described in the table
  	rtlsdr_i2c_write_reg(dev, i2c_addr, 0x01, i2c_value);
  	// Close the i2c_repeater
  	rtlsdr_set_i2c_repeater(dev, i2c_repeater_off);
  	// Reset the buffer
  	rtlsdr_reset_buffer(dev);
	}
}

void noise_collection(void)
{
	// Create thread for each RTL-SDRs data collection
	for(i = 0; i < NUM_SDRS; ++i)
	{
		if(pthread_create(&(sdrs[i].collection_t), NULL, collect_t, &i)) {
			fprintf(stderr, "Error creating thread\n");
			exit(1);
		}
	}

	// Wait until all threads have completed
	for(i = 0; i < NUM_SDRS; ++i)
	{
		if(pthread_join(sdrs[i].collection_t, NULL)) {
			fprintf(stderr, "Error joining thread\n");
			exit(1);
		}
	}
}

void rtlsdr_bias(void)
{
	for(i = 0; i < NUM_SDRS; ++i)
	{
		rtlsdr_dev_t *dev = sdrs[i].dev;

		// Set the bias tee by setting the gpio bit 0 to bias_off
  	rtlsdr_set_bias_tee(dev, bias_data);
  	// Set rtlsdr repeater for the i2communication via RTL2838
  	rtlsdr_set_i2c_repeater(dev, i2c_repeater_on);
  	// Set register to the output
  	rtlsdr_i2c_write_reg(dev, i2c_addr, 0x03, 00);
  	// Set value to the register as described in the table
  	rtlsdr_i2c_write_reg(dev, i2c_addr, 0x01, i2c_value);
  	// Close the i2c_repeater
  	rtlsdr_set_i2c_repeater(dev, i2c_repeater_off);
  	// Reset the buffer
  	rtlsdr_reset_buffer(dev);
	}
}

void data_collection(void)
{
	// Create thread for each RTL-SDRs data collection
	for(i = 0; i < NUM_SDRS; ++i)
	{
		if(pthread_create(&(sdrs[i].collection_t), NULL, collect_t, &i)) {
			fprintf(stderr, "Error creating thread\n");
			exit(1);
		}
	}

	// Wait until all threads have completed
	for(i = 0; i < NUM_SDRS; ++i)
	{
		if(pthread_join(sdrs[i].collection_t, NULL)) {
			fprintf(stderr, "Error joining thread\n");
			exit(1);
		}
	}
}
