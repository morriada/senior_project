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
#include "data_processing.h"

// Include external functions
extern void sdrs_setup(void);
extern void rtlsdr_setup(int, rtlsdr_dev_t *);
extern void rtlsdr_bias(int, uint8_t);
extern void collect(int, /*int,*/ rtlsdr_dev_t *);

// Initialize variables
extern pthread_mutex_t file;

int sdr0[2];
int sdr1[2];
int sdr2[2];

int m_time = 100000;

// DSP Thread
void * dodsp(void * ptr)
{
  // Initialize
  uint8_t * sdr0Data = malloc(2*BSIZE);
  uint8_t * sdr1Data = malloc(2*BSIZE);
  uint8_t * sdr2Data = malloc(2*BSIZE);

  memcpy(sdr0Data, sdrs[0].buffer[0], BSIZE);
  memcpy(sdr0Data + BSIZE, sdrs[0].buffer[1], BSIZE);
  memcpy(sdr1Data, sdrs[1].buffer[0], BSIZE);
  memcpy(sdr1Data + BSIZE, sdrs[1].buffer[1], BSIZE);
  memcpy(sdr2Data, sdrs[2].buffer[0], BSIZE);
  memcpy(sdr2Data + BSIZE, sdrs[2].buffer[1], BSIZE);

  fft_init();
  correlateInit();
  // Find Phase Difference
  findPhaseDifference(sdr0Data, sdr1Data, sdr2Data);
  // Find Signal in the Data Haystack
  DSP(sdr0Data, sdr1Data, sdr2Data);

  free(sdr0Data);
  free(sdr1Data);
  free(sdr2Data);

  pthread_exit(NULL);
}

// Individual SDR Collection Thread
void * collect_t(void * ptr)
{
  struct thread_struct * ts = (struct thread_struct *)ptr;

  // Tell Super thread we're at collection
  int val = 1;
  if (ts->id == 0) {
    write(sdr0[WRITE], &val, 1);
  } else if (ts->id == 1) {
    write(sdr1[WRITE], &val, 1);
  } else if (ts->id == 2) {
    write(sdr2[WRITE], &val, 1);
  }
  printf("here %d\n", ts->id);
  // Wait for Super thread to continue
  int ret = 1;
  while(ret)
  {
    if (ts->id == 0) {
      read(sdr0[READ], &ret, 1);
    } else if (ts->id == 1) {
      read(sdr1[READ], &ret, 1);
    } else if (ts->id == 2) {
      read(sdr2[READ], &ret, 1);
    }
  }
  // Collect Data
  collect(ts->id, /*ts->freq,*/ ts->dev);

  pthread_exit(NULL);
}

void * init_t(void * ptr)
{
  struct thread_struct * ts = (struct thread_struct *)ptr;
  rtlsdr_dev_t *dev = NULL;

  rtlsdr_open(&dev, ts->id);
  ts->dev = dev;
  sdrs[ts->id].dev = dev;
  rtlsdr_setup(ts->freq, dev);
  rtlsdr_reset_buffer(dev);
  rtlsdr_set_bias_tee(dev, 1);

  pthread_exit(NULL);
}


// Main Program
int main(void)
{
  // Declare variables
  int i, n;
  pthread_t dsp = (pthread_t)malloc(sizeof(pthread_t));
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
  n = 3;
  while(1)
  {
//    for(n = 0; n < 4; ++n)
//    {
      struct thread_struct tmp[3];

      for(i = 0; i < NUM_SDRS; ++i)
      {
        rtlsdr_dev_t *dev = NULL;
        tmp[i].id = i;
        tmp[i].freq = n;
        tmp[i].dev = dev;
        struct thread_struct * ts = &tmp[i];
        if(pthread_create(&(sdrs[i].initialize_t), NULL, init_t, (void *)ts)) {
          exit(1);
        }
      }

      for(i = 0; i < NUM_SDRS; ++i)
      {
        if(pthread_join(sdrs[i].initialize_t, NULL)) {
	         exit(1);
        }
      }

      rtlsdr_open(&(super.dev), 3);
      rtlsdr_bias(0, 0x1f);

      printf("Frequency: %d\n", freq[n]);

      // Create a collection thread for each RTL-SDR
      for(i = 0; i < NUM_SDRS; ++i)
      {
        struct thread_struct * ts = &tmp[i];
        if(pthread_create(&(sdrs[i].collection_t), NULL, collect_t, (void *)ts)) {
          fprintf(stderr, "Error creating thread\n");
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
      // Tell threads to continue
      ret = 0;
      write(sdr0[WRITE], &ret, 1);
      write(sdr1[WRITE], &ret, 1);
      write(sdr2[WRITE], &ret, 1);

      // Sleep for 100 milliseconds
      usleep(100000);
      // Switch RTL-SDRs Bias for Data Collection
      rtlsdr_bias(0, 0x00);

      // Wait for each collection thread to join
      for(i = 0; i < NUM_SDRS; ++i)
      {
        if(pthread_join(sdrs[i].collection_t, NULL)) {
          fprintf(stderr, "Error joining thread\n");
          exit(1);
        }
        // Close RTL-SDR device
        rtlsdr_close(sdrs[i].dev);
      }
      // Close Supervisory Channel
      rtlsdr_close(super.dev);

      // Perform DSP
      if(pthread_create(&dsp, NULL, dodsp, (void *)NULL)) {
        fprintf(stderr, "Error creating thread\n");
        exit(1);
      }
      if(pthread_join(dsp, NULL)) {
        fprintf(stderr, "Error joining thread\n");
        exit(1);
      }
 //   }
  }

  // Destroy mutex
  pthread_mutex_destroy(&file);

  // Free Space
  free_controls();

  return 0;
}
