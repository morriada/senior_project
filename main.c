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
extern void rtlsdr_setup(int);
extern void rtlsdr_freq(int, int);
extern void rtlsdr_bias(uint8_t);
extern void data_collection(int);

// DSP Thread
void dodsp(void * ptr)
{
  //
}

// Individual SDR Collection Thread
void * collect_t(void * ptr)
{
  struct thread_struct * ts = (struct thread_struct *)ptr;

  // Open RTL-SDR device
  rtlsdr_open(&(sdrs[ts->id].dev), sdrs[ts->id].id);
  // Setup RTL-SDRs
  rtlsdr_setup(ts->id);
  // Set RTL-SDRs for desired frequency
  rtlsdr_freq(ts->id, ts->freq);
  // Collect Data
  data_collection(ts->id);
  // Close RTL-SDR device
  rtlsdr_close(sdrs[ts->id].dev);

  return NULL;
}

// Supervisory SDR Thread
void * super_t(void * ptr)
{
  int i;
  int * n = (int *)ptr;

  // Set RTL-SDRs Bias for Noise Collection
  rtlsdr_bias(0x1f);

  // Create a collection thread for each RTL-SDR  
  for(i = 0; i < NUM_SDRS; ++i)
  {
    struct thread_struct * ts;
    ts->id = i;
    ts->freq = *n;
    if(pthread_create(&(sdrs[i].collection_t), NULL, collect_t, &ts)) {
      fprintf(stderr, "Error creating thread\n");
      exit(1);
    }
  }
  
  // Sleep for 100 milliseconds
  sleep(0.1);
  // Switch RTL-SDRs Bias for Data Collection
  rtlsdr_bias(0x00);

  // Wait for each collection thread to join
  for(i = 0; i < NUM_SDRS; ++i)
  {
    if(pthread_join(sdrs[i].collection_t, NULL)) {
      fprintf(stderr, "Error joining thread\n");
      exit(1);
    }
  }

  // Create a DSP thread
  
  return NULL;
}


// Main Program
int main(void)
{
  // Declare variables
  int n;
  // Prepare structures
  sdrs_setup();

  while(1)
  {
    for(n = 0; n < 4; ++n)
    {
      if(pthread_create(&(super.collection_t), NULL, super_t, &n)) {
        fprintf(stderr, "Error creating thread\n");
        exit(1);
      }

      if(pthread_join(super.collection_t, NULL)) {
        fprintf(stderr, "Error joining thread\n");
        exit(1);
      }
    }
  }

  return 0;
}
