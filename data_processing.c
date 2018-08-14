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
 * GNU General Public Licens e for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* Standard Includes */
#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>

/* Signal Processing Includes */
#include "fftw3.h"

static int nreceivers = 0, corrlen = 0, fft1n = 0, fft2n = 0;

static int sync_debug = 1;

static fftwf_complex *fft1in, *fft1out, *fft2in, *fft2out;
static fftwf_plan fft1plan, fft2plan;

int fft_init(int numReceivers, int sync_len){
	int i, arraysize;
	nreceivers = numReceivers;
	corrlen = sync_len;
	if(nreceivers > NRECEIVERS_MAX) return -1;


	/* Half of FFT input will be zero padded
	   this is for cross correlation		*/
	fft1n = corrlen*2;
	fft2n = fft1n * CORRELATION_OVERSAMPLE;

	arraysize = nreceivers * fft1n;
	fft1in = fftwf_mallloc(arraysize * sizeof(*fft1in));
	fft1out = fftwf_mallloc(arraysize * sizeof(*fft1out));
	for(i = 0; i < arraysize; i++){
		fft1in[i] = fft1out[i] = 0;
	}

	arraysize = (nreceiver - 1) *fft2n;
	fft2in = fftwf_mallloc(arraysize * sizeof(*fft2in));
	fft2out = fftwf_mallloc(arraysize * sizeof(*fft2out));
	for(i = 0; i < arraysize; i++){
		fft2in[i] = fft1out[i] = 0;
	}

	fft1plan = fftwf_plan_many_dft(
		1, &fft1n, nreceiver,
		fft1in, NULL, 1, fft1n,
		fft1out, NULL, 1, fft1n,
		FFTW_FORWARD, FFTW_ESTIMATE
	);

	fft2plan = fftwf_plan_many_dft(
		1, &fft2n, nreceiver,
		fft2in, NULL, 1, fft2n,
		fft2out, NULL, 1, fft2n,
		FFTW_FORWARD, FFTW_ESTIMATE
	);

	return 0;
}

int sync_block(int blocksize, csample_t **buffers, float *timediffs, float *phasediffs) {
	int ri, i;
	if(blocksize < corrlen) return -1;
	
	for(ri = 0; ri < nreceivers; ri++) {
		csample_t *buf = buffers[ri];
		fftwf_complex *floatbuf = fft1in + ri * fft1n;
		for(i = 0; i < corrlen; i++)
			floatbuf[i] = (buf[i][0] + I*buf[i][1]) - (127.4f+127.4f*I);
	}
	fftwf_execute(fft1plan);

	for(ri = 1; ri < nreceivers; ri++) {
		/* cross correlation of receiver number ri and receiver 0 */
		fftwf_complex *f1o = fft1out + ri * fft1n;
		fftwf_complex *f2i = fft2in  + (ri-1) * fft2n;
		/* positive frequencies: */
		for(i = 0; i < fft1n/2; i++)
			f2i[i] = fft1out[i] * conjf(f1o[i]);
		/* negative frequencies: */
		f2i += fft2n - fft1n;
		for(i = fft1n/2; i < fft1n; i++)
			f2i[i] = fft1out[i] * conjf(f1o[i]);
	}
	fftwf_execute(fft2plan);
	
	timediffs[0] = 0;
	phasediffs[0] = 0;
	if(sync_debug) 
		printf("%d ",(int)time(NULL));
	for(ri = 1; ri < nreceivers; ri++) {
		fftwf_complex *f2o = fft2out  + (ri-1) * fft2n;
		float maxmagsq = 0, phasedifference;
		float timedifference = 0;
		float complex maxc = 0;
		float y1, y2, y3;
		int maxi = 1;
		for(i = 0; i < fft2n; i++) {
			float complex c = f2o[i];
			float magsq = crealf(c)*crealf(c) + cimagf(c)*cimagf(c);
			if(magsq > maxmagsq) {
				maxmagsq = magsq;
				maxc = c;
				maxi = i;
			}
		}
		/* parabolic interpolation for more fractional sample resolution
		(math from http://dspguru.com/dsp/howtos/how-to-interpolate-fft-peak ) */
		y1 = cabsf(f2o[(maxi-1 + fft2n) % fft2n]);
		y2 = cabsf(maxc);
		y3 = cabsf(f2o[(maxi+1 + fft2n) % fft2n]);
		//printf("%E %E %E   ", y1, y2, y3);

		if(maxi >= fft2n/2) maxi -= fft2n;
		timedifference = maxi;
		timedifference += (y3-y1) / (2*(2*y2-y1-y3));
		timedifference *= 1.0 / CORRELATION_OVERSAMPLE;
		timediffs[ri] = timedifference;
		
		phasediffs[ri] = phasedifference = cargf(maxc);
		if(sync_debug)
			printf("%9.2f %E %6.2f  ", timedifference, maxmagsq, 57.2957795 * phasedifference);
	}
	if(sync_debug) {
		printf("\n");
		fflush(stdout);
	}
	return 0;
}


/*Correlate Data */
uint8_t correlateData(){
	//Take Cross Correlate Input Data
}

/*FFT Data*/
uint8_t * fftData(fftwf_complex *fftin){
	//Take FFTs of Data
	fftplan = fftwf_plan_many_dft(1, &fftplan, )
}

/*Find Peaks*/
uint8_t * findpeaks(){



}