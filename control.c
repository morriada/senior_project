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
volatile uint32_t freq[] = {
	148000000,
	149000000,
	150000000,
	151000000
};

int i, c_time = 50000;
int auto_gain = 0;
int disable_dither = 0;
int i2c_repeater_on = 1;
int i2c_repeater_off = 0;
uint32_t i2c_addr = 0x40;
uint32_t sample_rate = 2400000;
uint32_t if_freq = 7200000;

void sdrs_setup(void)
{
	for(i = 0; i < NUM_SDRS; ++i)
	{
		// Initialize Variables
		sdrs[i].id = i;
		sdrs[i].blocksize = BSIZE;
		sdrs[i].buffer[0] = calloc(BSIZE, sizeof(uint8_t));
		sdrs[i].buffer[1] = calloc(BSIZE, sizeof(uint8_t));
		sdrs[i].collection_t = (pthread_t)malloc(sizeof(pthread_t));
		sdrs[i].initialize_t = (pthread_t)malloc(sizeof(pthread_t));
	}
	super.id = NUM_SDRS;
	super.blocksize = 0;
}

void free_controls(void)
{
	for(i = 0; i < NUM_SDRS; ++i)
	{
		free((void *)sdrs[i].buffer[0]);
		free((void *)sdrs[i].buffer[1]);
		free((void *)sdrs[i].collection_t);
		free((void *)sdrs[i].initialize_t);
	}
}

void rtlsdr_setup(int f, rtlsdr_dev_t *dev)
{
	int r;
	// Set the sample rate of the rtl-sdr
	if((r = rtlsdr_set_sample_rate(dev, sample_rate)) < 0)
		printf("WARNING: [%d] Failed to set sample rate.\n", r);
	// Disable dithering
	if((r = rtlsdr_set_dithering(dev, disable_dither)) < 0)
		printf("WARNING: [%d] Failed to set dithering.\n", r);
	// Set the IF frequency
//	if((r = rtlsdr_set_if_freq(dev, if_freq)) < 0)
//		printf("WARNING: [%d] Failed to set if frequency.\n", r);
	// Set the center frequency
	if((r = rtlsdr_set_center_freq(dev, freq[f])) < 0)
		printf("WARNING: [%d] Failed to set if frequency.\n", r);
	// Set the tuner gain mode to automatic
	if((r = rtlsdr_set_tuner_gain_mode(dev, auto_gain)) < 0)
		printf("WARNING: [%d] Failed to set tuner gain.\n", r);
}

void file_save(int frequency, int band, float phase, int mort)
{
	FILE *fp;

	// Current Time
	char path[SIZE], line[SIZE];
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

	// Save buffer to file for sdr_num and frequency
	snprintf(path, sizeof(char) * SIZE, "/home/pi/data/freq%i.csv",
					freq[frequency] + (band * 100000));
	snprintf(line, sizeof(char) * SIZE, "%i %l %i\n",
					freq[frequency] + (band * 100000), phase, mort);
	fp = fopen(path, "wb");
	fwrite(line, 1, BSIZE, fp);
	fclose(fp);
}

/*void file_save(int sdr_num, int f)
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

	// Save buffer to file for sdr_num and frequency
	snprintf(path, sizeof(char) * SIZE, "/home/pi/data/sdr%i_freq%i_%i%i%i%i%i%i.dat",
					sdr_num, freq[f], yr, mon, day, hr, min, sec);
	pthread_mutex_lock(&file);
	fp = fopen(path, "wb");
	fwrite(sdrs[sdr_num].buffer[0], 1, BSIZE, fp);
	fwrite(sdrs[sdr_num].buffer[1], 1, BSIZE, fp);
	fclose(fp);
	pthread_mutex_unlock(&file);
}*/

void collect(int id, int f, rtlsdr_dev_t *dev)
{
	// Initialize collection variables
	int ret, blocksize, n_read, idx;

	// Collect data from all RTL-SDRs
	for(idx = 0; idx < 2; ++idx){
		ret = n_read = 0;
		blocksize = sdrs[id].blocksize;

		ret = rtlsdr_read_sync(dev, sdrs[id].buffer[idx], blocksize, &n_read);

		// Check for errors
		if(ret < 0) {
			fprintf(stderr, "Runtime error: %d at %s:%d\n", ret, __FILE__, __LINE__);
		} else if(n_read < blocksize) {
			fprintf(stderr, "Short read sdr: %d: %d/%d ret: %d\n", id, n_read, blocksize, ret);
		} else {
			fprintf(stderr, "Read %d:%d\n", id, idx);
		}
	}

	// Save data to file
	//file_save(id, f);
}

void rtlsdr_bias(int bias, uint8_t i2c_val)
{
	int r;
	// Set the bias tee by setting the gpio bit 0 to bias_off
	if((r = rtlsdr_set_bias_tee(super.dev, bias)) < 0)
		printf("WARNING: [%d] Failed to set bias tee.\n", r);
	// Set rtlsdr repeater for the i2communication via RTL2838
	rtlsdr_set_i2c_repeater(super.dev, i2c_repeater_on);
	// Set register to the output
	if((r = rtlsdr_i2c_write_reg(super.dev, i2c_addr, 0x03, 0x00)) < 0)
		printf("WARNING: [%d] Failed to write to i2c.\n", r);
	// Set value to the register as described in the table
	if((r = rtlsdr_i2c_write_reg(super.dev, i2c_addr, 0x01, i2c_val)) < 0)
		printf("WARNING: [%d] Failed to write to i2c.\n", r);
	// Close the i2c_repeater
	rtlsdr_set_i2c_repeater(super.dev, i2c_repeater_off);
}
