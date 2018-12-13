/*
 * Data Processing - details
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

#ifndef SENIOR_PROJECT_DATA_H
#define SENIOR_PROJECT_DATA_H

/* Standard Includes */
#include <stdint.h>
#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* Signal Processing Includes */
#include "fftw3.h"

#define CALIBRATION_LENGTH 660000
#define CALIBRATION_SKIP 400000
#define FFT_TIME 5
#define SAMPLE_TIME 2130
#define SAMPLE_LENGTH 5279744
#define SAMPLE_FREQUENCY 1200000
#define SAMPLE_RATE 2400000
#define COLLAR_OFFSET 43000
#define COLLAR_TOLERANCE 2500
#define NUM_BANDS 12


float angleOfArrival[NUM_BANDS];

/*
 * correlateInit - Description.
 */
void correlateInit(void);

/*
 * fft_init - Description.
 */
void fft_init(void);

/*
 * DSP - Description.
 * @param SDR1_data Data from first SDR
 * @param SDR2_data Data from second SDR
 * @param SDR3_data Data from third SDR
 */
int DSP(uint8_t *SDR1_data, uint8_t *SDR2_data, uint8_t *SDR3_data);

/*
 * findPhaseDifference - Description.
 * @param parameter Description
 */
void findPhaseDifference(uint8_t *SDR1_cal, uint8_t *SDR2_cal, uint8_t *SDR3_cal);

#endif // SENIOR_PROJECT_CONTROL_H
