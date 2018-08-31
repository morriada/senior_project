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
extern void sdrs_setup(void);
extern void rtlsdr_setup(int, int);
extern void rtlsdr_bias(int, uint8_t);
extern void collect(int, int);

// Initialize variables
extern pthread_mutex_t file;

int sdr0[2];
int sdr1[2];
int sdr2[2];

int m_time = 1000000;

// DSP Thread
void dodsp(void * ptr)
{
  //
}

// Individual SDR Collection Thread
void * collect_t(void * ptr)
{
  struct thread_struct * ts = (struct thread_struct *)ptr;
  int r;

  rtlsdr_open(&(sdrs[ts->id].dev), ts->id);
  rtlsdr_setup(ts->id, ts->freq);
  rtlsdr_reset_buffer(sdrs[ts->id].dev);
  rtlsdr_set_bias_tee(sdrs[ts->id].dev, 1);

  // Tell Super thread we're at collection
  int val = 1;
  if (ts->id == 0) {
    write(sdr0[WRITE], &val, 1);
  } else if (ts->id == 1) {
    write(sdr1[WRITE], &val, 1);
  } else if (ts->id == 2) {
    write(sdr2[WRITE], &val, 1);
  }

  // Wait for Super thread to continue
  int ret = 0;
  if (ts->id == 0) {
    read(sdr0[WRITE], &ret, 1);
  } else if (ts->id == 1) {
    read(sdr1[WRITE], &ret, 1);
  } else if (ts->id == 2) {
    read(sdr2[WRITE], &ret, 1);
  }

  // Collect Data
  collect(ts->id, ts->freq);

  // Close RTL-SDR device
  rtlsdr_close(sdrs[i].dev);

  pthread_exit(NULL);
}


// Main Program
int main(void)
{
  // Declare variables
  int i, n, r;
  // Prepare structures
  sdrs_setup();
  // Initialize pipes
  if((pipe(sdr0) < 0) || (pipe(sdr1) < 0) || (pipe(sdr2) < 0))
  {
    printf("\n pipe init has failed\n");
    return 1;
  }
  // Initialize mutex
  if(pthread_mutex_init(&file, NULL))
   {
    printf("\n mutex init has failed\n");
    return 1;
  }

  while(1)
  {
    for(n = 0; n < 4; ++n)
    {
      struct thread_struct tmp[3];

      // Create a collection thread for each RTL-SDR
      for(i = 0; i < NUM_SDRS; ++i)
      {
        tmp[i].id = i;
        tmp[i].freq = n;
        struct thread_struct * ts = &tmp[i];
        if(pthread_create(&(sdrs[i].collection_t), NULL, collect_t, (void *)ts)) {
          //fprintf(stderr, "Error creating thread\n");
          exit(1);
        }
      }

      // Wait for SDRs to be at collection
      int ret;
      for(i = 0; i < NUM_SDRS; ++i)
      {
        ret = 0;
        while(!ret)
        {
          if (i == 0) {
            read(sdr0[READ], &ret, 1);
          } else if (i == 1) {
            read(sdr1[READ], &ret, 1);
          } else if (i == 2) {
            read(sdr2[READ], &ret, 1);
          }
        }
      }

      // Open Supervisory SDR
      rtlsdr_open(&(super.dev), 3);
      // Change Bias Tee
      rtlsdr_bias(0, 0x1f);

      // Tell threads to continue
      ret = 1;
      for(i = 0; i < NUM_SDRS; ++i)
      {
        if (i == 0) {
          write(sdr0[WRITE], &ret, 1);
        } else if (i == 1) {
          write(sdr1[WRITE], &ret, 1);
        } else if (i == 2) {
          write(sdr2[WRITE], &ret, 1);
        }
      }

      // Sleep for 100 milliseconds
      usleep(m_time);
      // Switch RTL-SDRs Bias for Data Collection
      rtlsdr_bias(0, 0x00);

      // Wait for each collection thread to join
      for(i = 0; i < NUM_SDRS; ++i)
      {
        if(pthread_join(sdrs[i].collection_t, NULL)) {
          //fprintf(stderr, "Error joining thread\n");
          exit(1);
        }
      }
      // Close Supervisory Channel
      rtlsdr_close(super.dev);
    }
  }

  // Destroy mutex
  pthread_mutex_destroy(&file);

  return 0;
}
