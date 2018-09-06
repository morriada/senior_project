/*
 * control.h - Initialize peripherals and contol their operations.
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

#ifndef SENIOR_PROJECT_CONTROL_H
#define SENIOR_PROJECT_CONTROL_H

// Include Standard Libraries
#include <errno.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

// Include SDR Library
#include "rtl-sdr.h"

// Define Macros
//#define BSIZE 10559488
//#define BSIZE 499712
#define BSIZE 5279744
#define NUM_SDRS 3
#define READ 0
#define SIZE 100
#define WRITE 1

// Declare variables and structures
static pthread_mutex_t file;

typedef struct rtlsdr_struct {
	int blocksize;
	uint32_t id;
	uint8_t *buffer[2];
	rtlsdr_dev_t *dev;
	pthread_t collection_t, initialize_t;
} rtlsdr_struct;
struct rtlsdr_struct super;
struct rtlsdr_struct sdrs[3];

typedef struct thread_struct {
	int id;
	int freq;
	rtlsdr_dev_t *dev;
} thread_struct;

/*
 * sdrs_setup - Initializes the structures for each RTL-SDR.
 */
void sdrs_setup(void);

/*
 * rtlsdr_setup - Sets up an individual RTL-SDR at the beginning of the program.
 * @param id ID of RTL-SDR - expecting an integer from 0 to NUM_SDRS
 */
void rtlsdr_setup(int id, int f, rtlsdr_dev_t *dev);

/*
 * rtlsdr_bias - Sets bias of Supervisory RTL-SDR for data collection from the
 * 							noise card or from the antenna's.
 * @param i2c_val 2 byte value for register on RTL-SDR
 */
void rtlsdr_bias(int bias, uint8_t i2c_val);

/*
 * collect - Colects the data as set up from the rtlsdr_setup function.
 * @param id ID of RTL-SDR - expecting an integer from 0 to NUM_SDRS
 * @param f frequency id - expecting an integer from 0 to 3
 */
void collect(int id, int f, rtlsdr_dev_t *dev);

#endif // SENIOR_PROJECT_CONTROL_H
