/*
 * control.c - Initialize peripherals and contol their operations.
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

int i;
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

void sdrs_setup(void )
{
	for(i = 0; i < NUM_SDRS; ++i)
	{
		// Initialize Variables
		sdrs[i].id = i;
		sdrs[i].blocksize = BSIZE;
		sdrs[i].calibration = 1;
	}
}

void rtlsdr_setup(rtlsdr_struct *sdr, int f)
{
	// Set the sample rate of the rtl-sdr
	rtlsdr_set_sample_rate(sdr->dev, sample_rate);
	// Disable dithering
	rtlsdr_set_dithering(sdr->dev, disable_dither);
	// Set the tuner gain mode to automatic
	rtlsdr_set_tuner_gain_mode(sdr->dev, auto_gain);
}

void rtlsdr_freq(rtlsdr_struct *sdr, int f)
{
	// Set the IF frequency
	rtlsdr_set_if_freq(sdr->dev, if_freq);
	// Set the center frequency
	rtlsdr_set_center_freq(sdr->dev, freq[f]);
}

void file_save(int sdr_num)
{
	FILE *fp;
	
	// Current Time
	char path[SIZE];
	time_t curtime;
	struct tm *loctime;
	curtime = time(NULL);
	loctime = localtime(&curtime);
	int sec = loctime->tm_sec;
	int min = loctime->tm_min;
  	int hr = loctime->tm_hour;
  	int day = loctime->tm_mday;
  	int mon = loctime->tm_mon;
  	int yr = loctime->tm_year;

	// Save buffer to file for sdr_num
	if(sdrs[sdr_num].calibration) {
		sdrs[sdr_num].calibration = 0;
		// Save buffer to Calibration file
		snprintf(path, sizeof(char) * SIZE, "/home/pi/data/calibration%i_%i%i%i%i%i%i.dat",
						sdr_num, yr, mon, day, hr, min, sec);
		fp = fopen(path, "w");
		int z;
		// Write to file byte by byte
		for(z = 0; z < BSIZE; ++z)
			fwrite(&sdrs[sdr_num].buffer[z], 1, 1, fp);
		fclose(fp);
	} else {
		sdrs[sdr_num].calibration = 1;
		// Save buffer to Data file
		snprintf(path, sizeof(char) * SIZE, "/home/pi/data/data%i_%i%i%i%i%i%i.dat",
						sdr_num, yr, mon, day, hr, min, sec);
		fp = fopen(path, "w");
		int z;
		// Write to file byte by byte
		for(z = 0; z < BSIZE; ++z)
			fwrite(&sdrs[sdr_num].buffer[z], 1, 1, fp);
		fclose(fp);
	}
}

void collect(rtlsdr_struct *sdr, int f)
{
	// Initialize collection variables
	int ret, blocksize, n_read;

	// Collect data from all RTL-SDRs
	ret = n_read = 0;
	blocksize = sdr->blocksize;
	ret = rtlsdr_read_sync(sdr->dev, sdr->buffer, blocksize, &n_read);

	// Check for errors
	if(ret < 0) {
		fprintf(stderr, "Runtime error: %d at %s:%d\n", ret, __FILE__, __LINE__);
	} else if(n_read < blocksize) {
		fprintf(stderr, "Short read %d: %d/%d\n", sdr->id, n_read, blocksize);
	} else {
		fprintf(stderr, "Read %d\n", sdr->id);
	}

	// Save data to file
	file_save(i);
}

void rtlsdr_calibration(rtlsdr_struct *sdr, int f)
{
	// Set the bias tee by setting the gpio bit 0 to bias_off
  	rtlsdr_set_bias_tee(sdr->dev, bias_noise);
  	// Set rtlsdr repeater for the i2communication via RTL2838
  	//rtlsdr_set_i2c_repeater(dev, i2c_repeater_on);
  	// Set register to the output
  	//rtlsdr_i2c_write_reg(dev, i2c_addr, 0x03, 00);
  	// Set value to the register as described in the table
  	//rtlsdr_i2c_write_reg(dev, i2c_addr, 0x01, i2c_value);
  	// Close the i2c_repeater
  	//rtlsdr_set_i2c_repeater(dev, i2c_repeater_off);
  	// Reset the buffer
  	rtlsdr_reset_buffer(sdr->dev);
}

void noise_collection(rtlsdr_struct *sdr, int f)
{
	collect(sdr, f);
	sdr->calibration = 0;
}

void rtlsdr_bias(rtlsdr_struct *sdr, int f)
{
	// Set the bias tee by setting the gpio bit 0 to bias_off
  	rtlsdr_set_bias_tee(sdr->dev, bias_data);
  	// Set rtlsdr repeater for the i2communication via RTL2838
  	//rtlsdr_set_i2c_repeater(dev, i2c_repeater_on);
  	// Set register to the output
  	//rtlsdr_i2c_write_reg(dev, i2c_addr, 0x03, 00);
  	// Set value to the register as described in the table
  	//rtlsdr_i2c_write_reg(dev, i2c_addr, 0x01, i2c_value);
  	// Close the i2c_repeater
  	//rtlsdr_set_i2c_repeater(dev, i2c_repeater_off);
  	// Reset the buffer
  	rtlsdr_reset_buffer(sdr->dev);
}

void data_collection(rtlsdr_struct *sdr, int f)
{
	collect(sdr, f);
	sdr->calibration = 1;
}
