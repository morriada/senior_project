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
#define NUM_SDRS 3
#define SIZE 100

// Declare variables and structures
typedef struct rtlsdr_struct {
	int blocksize, calibration;
	uint32_t id;
	uint8_t *buffer;
	rtlsdr_dev_t *dev;
	pthread_t collection_t;
} rtlsdr_struct;
struct rtlsdr_struct sdrs[3];

/*
 * usage - Prints usage description of program
 */
void usage(void);

/*
 * rtlsdr_setup - Sets up an individual RTL-SDR at the beginning of
 *                the program.
 */
void rtlsdr_setup(void);

/*
 *
 */
void rtlsdr_calibration(void);

/*
 *
 */
void noise_collection(void);

/*
 *
 */
void rtlsdr_bias(void);

/*
 *
 */
void data_collection(void);

#endif // SENIOR_PROJECT_CONTROL_H
