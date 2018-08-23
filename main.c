/*
 * Tracking Tower - Recieves and Processes data from SDRs.
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

// Include Standard Libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Include SDR Library
#include "rtl-sdr.h"

// Include Project Libraries
#include "control.h"
//#include "data_processing.h"

// Include external functions
extern void usage(void);
extern void sdrs_setup(void);
extern void rtlsdr_setup(rtlsdr_struct *, int);
extern void rtlsdr_freq(rtlsdr_struct *, int);
extern void noise_collection(rtlsdr_struct *, int);
extern void rtlsdr_bias(rtlsdr_struct *, int);
extern void data_collection(rtlsdr_struct *, int);

int main(/*int argc, char** argv*/)
{
  // Declare variables
  int i, n;
  // Prepare structures
  sdrs_setup();

  while(1)
  {
    for(n = 0; n < 4; ++n)
    {
      for(i = 0; i < NUM_SDRS; ++i)
      {
        // Open RTL-SDR device
        rtlsdr_open(&(sdrs[i].dev), sdrs[i].id);
        // Setup RTL-SDRs
        rtlsdr_setup(&sdrs[i], n);
        // Set RTL-SDRs for desired frequency
        rtlsdr_freq(&sdrs[i], n);
        // Prepare RTL-SDRs for Calibration
        rtlsdr_calibration(&sdrs[i], n);
        // Collect Data for Calibration
        noise_collection(&sdrs[i], n);
        // Switch RTL-SDRs Bias for Data Collection
        rtlsdr_bias(&sdrs[i], n);
        // Collect Accurate Data
        data_collection(&sdrs[i], n);
        // Close RTL-SDR device
        rtlsdr_close(sdrs[i].dev);
      }
    }
  }

  return 0;
}
