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
extern void rtlsdr_setup(void);
extern void rtlsdr_freq(int);
extern void noise_collection(void);
extern void rtlsdr_bias(void);
extern void data_collection(void);

int main(/*int argc, char** argv*/)
{
  // Declare variables
  int n;
  // Setup RTL-SDRs
  rtlsdr_setup();

  while(1)
  {
    for(n = 0; n < 4; ++n)
    {
      // Set RTL-SDRs for desired frequency
      rtlsdr_freq(n);
      // Prepare RTL-SDRs for Calibration
      rtlsdr_calibration();
      // Collect Data for Calibration
      noise_collection();
      // Switch RTL-SDRs Bias for Data Collection
      rtlsdr_bias();
      // Collect Accurate Data
      data_collection();
      // Find phase difference from Calibration Data
      //function
      // Find Peaks within each band
      //function
      // Correct for Phase Difference
      //function
      // Determine Direction of Arrival
      //  with Phase Interferometry
      //function
      // Save data to files
      //function
    }
  }

  return 0;
}
