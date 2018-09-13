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

#include "data_processing.h"

#define numFFTs ( SAMPLE_TIME / FFT_TIME )
#define distance ( SAMPLE_LENGTH / numFFTs )

static int nreceivers = 0, corrlen = 0, fft1n = 0, fft2n = 0;

static int sync_debug = 1;

static fftw_complex *fft1in, *fft1out, *fft2in, *fft2out, *fft3in, *fft3out;
static fftw_plan fft1plan, fft2plan, fft3plan;
float complex peaks[3][NUM_BANDS][numFFTs];
float angleOfArrival[10], phaseCorrectionAC, phaseCorrectionBC;

/*int corr_init(int numReceivers, int sync_len){
	int i, arraysize;
	nreceivers = numReceivers;
	corrlen = sync_len;
	if(nreceivers > NRECEIVERS_MAX) return -1;


	// Half of FFT input will be zero padded
	//   this is for cross correlation
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
}*/


/*void dsp(fftw_complex *in, fftw_complex *out, int N){
	fftw_plan p;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
}*/

/* Initialize Correlation */
void correlateInit(void)
{
    phaseCorrectionAC = 0.0;
    phaseCorrectionBC = 0.0;
}

/*Correlate Data */
void xcorr(fftw_complex * signala, fftw_complex * signalb, fftw_complex * result, int N)
{
    fftw_complex * signala_ext = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (2 * N - 1));
    fftw_complex * signalb_ext = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (2 * N - 1));
    fftw_complex * out_shifted = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (2 * N - 1));
    fftw_complex * outa = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (2 * N - 1));
    fftw_complex * outb = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (2 * N - 1));
    fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (2 * N - 1));

    fftw_plan pa = fftw_plan_dft_1d(2 * N - 1, signala_ext, outa, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan pb = fftw_plan_dft_1d(2 * N - 1, signalb_ext, outb, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan px = fftw_plan_dft_1d(2 * N - 1, out, result, FFTW_BACKWARD, FFTW_ESTIMATE);

    //zeropadding
    memset (signala_ext, 0, sizeof(fftw_complex) * (N - 1));
    memcpy (signala_ext + (N - 1), signala, sizeof(fftw_complex) * N);
    memcpy (signalb_ext, signalb, sizeof(fftw_complex) * N);
    memset (signalb_ext + N, 0, sizeof(fftw_complex) * (N - 1));

    fftw_execute(pa);
    fftw_execute(pb);

    fftw_complex scale = 1.0/(2 * N -1);
    for (int i = 0; i < 2 * N - 1; i++)
        out[i] = outa[i] * conj(outb[i]) * scale;

    fftw_execute(px);

    fftw_destroy_plan(pa);
    fftw_destroy_plan(pb);
    fftw_destroy_plan(px);

    fftw_free(signala_ext);
    fftw_free(signalb_ext);
    fftw_free(out_shifted);
    fftw_free(out);
    fftw_free(outa);
    fftw_free(outb);

    fftw_cleanup();

    return;
}

/*FFT Data*/
/*uint8_t * fftData(fftwf_complex *fftin){
	//Take FFTs of Data
    fftplan = fftwf_plan_many_dft(1, &fftplan, );
}*/

/*Find Peaks*/
uint8_t * findpeaks(void)
{
    return 0;
}

float phaseInterferometry(int band, int location)
{
    int i;
    float phase[3], top, bottom;
    for(i = 0; i < 3; i++){
        phase[i] = (atanf(cimagf(peaks[i][band][location])/crealf(peaks[i][band][location]))) * (180.0/M_PI);
        if(cimagf(peaks[i][band][location]) > 0 && crealf(peaks[i][band][location]) < 0)
            phase[i] = (phase[i] * -1.0) + 90.0;
        else if(cimagf(peaks[i][band][location]) < 0 && crealf(peaks[i][band][location]) < 0)
            phase[i] = (phase[i] * -1.0) - 90.0;
    }
    top = 1.1547 * (phase[0] - phase[2] + phaseCorrectionAC);
    bottom = phase[1] - phase[2] + phaseCorrectionBC;

    if(top > 180.0)
        top -= 360.0;
    else if(top < -180.0)
        top += 360.0;
    if(bottom > 180.0)
        bottom -= 360.0;
    else if(bottom < -180.0)
        bottom += 360.0;

    return atanf((top/bottom) - 0.57735)*(180/M_PI);

}

float findSignal(void)
{
    //Initializ Variables
    int i, j, k, maxLocation[3][10] = {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};
    float max[3][10] = {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};
    float averageNoise[3][10] = {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};

    //Find the max of each band, location of the max, and the average noise level in that band.
    for(i = 0; i < 3; i++){
        for(j = 0; j < NUM_BANDS; j++){
            for(k = 0; k < numFFTs; k++){
                float peakmag = sqrt(crealf(peaks[i][j][k])*crealf(peaks[i][j][k]) + cimagf(peaks[i][j][k])*cimagf(peaks[i][j][k]));
                averageNoise[i][j] += peakmag;
                if(max[i][j] < peakmag){
                    max[i][j] = peakmag;
                    maxLocation[i][j] = k;
                }
            }
            averageNoise[i][j] /= (float)numFFTs;
        }
    }

    //for(i = 0; i< 10; i++)
    //    printf("Max Location of Band %d: %f\n", i, max[0][i]);
    int numLeft[3][10], numRight[3][10];
    int badLeftFlag = 0, badRightFlag = 0, stopLeft = 0, stopRight = 0;

    //Now we verify that the signal exists
    for(i = 0; i < 3; i++){
        for(j = 0; j < 10; j++){
            stopLeft = 0;
            stopRight = 0;
            numLeft[i][j] = 0;
            numRight[i][j] = 0;
            badRightFlag = 0;
            badLeftFlag = 0;
            float threshold = max[i][j] - ( (max[i][j]-averageNoise[i][j]) / 2 );
            //printf("SDR %d Band %d Threashold: %f\n",i,j, threshold);
            for(k = 1; k < 4; k++){
                float peakmagMin = sqrt(crealf(peaks[i][j][maxLocation[i][j] - k])*crealf(peaks[i][j][maxLocation[i][j] - k]) + cimagf(peaks[i][j][maxLocation[i][j] - k])*cimagf(peaks[i][j][maxLocation[i][j] - k]));
                float peakmagPlus = sqrt(crealf(peaks[i][j][maxLocation[i][j] + k])*crealf(peaks[i][j][maxLocation[i][j] + k]) + cimagf(peaks[i][j][maxLocation[i][j] + k])*cimagf(peaks[i][j][maxLocation[i][j] + k]));
                //if(j == 4)
                //    printf("%f\n",peakmagMin);
                if(peakmagMin > threshold && !stopLeft){
                    numLeft[i][j] = k;
                }

                else{
                    if(!badLeftFlag)
                        badLeftFlag = 1;
                    else
                        stopLeft = 1;
                }
                if(peakmagPlus > threshold && !stopRight){
                    numRight[i][j] = k;
                }
                else{
                    if(!badRightFlag)
                        badRightFlag = 1;
                    else
                        stopRight = 1;
                }
            }
        }
    }

    //Pick the points to use
    int locationForAnalysis[10];


    for(j = 0; j < 10; j++){
        if(((maxLocation[0][j] - numLeft[0][j]) >= (maxLocation[1][j] - numLeft[1][j])) && ((maxLocation[0][j] - numLeft[0][j]) >= (maxLocation[2][j] - numLeft[2][j])))
            locationForAnalysis[j] = maxLocation[0][j] - numLeft[0][j];
        else if(((maxLocation[1][j] - numLeft[1][j]) >= (maxLocation[0][j] - numLeft[0][j])) && ((maxLocation[1][j] - numLeft[1][j]) >= (maxLocation[2][j] - numLeft[2][j])))
            locationForAnalysis[j] = (maxLocation[1][j] - numLeft[1][j]);
        else if(((maxLocation[2][j] - numLeft[2][j]) >= (maxLocation[0][j] - numLeft[0][j])) && ((maxLocation[2][j] - numLeft[2][j]) >= (maxLocation[1][j] - numLeft[1][j])))
            locationForAnalysis[j] = (maxLocation[2][j] - numLeft[2][j]);

        if( ((numLeft[0][j] + numRight[0][j] + 1) >= 3) && ((numLeft[0][j] + numRight[0][j] + 1) <= 5) && ((numLeft[1][j] + numRight[1][j] + 1) >= 3) && ((numLeft[1][j] + numRight[1][j] + 1) <= 5) &&((numLeft[2][j] + numRight[2][j] + 1) >= 3) && ((numLeft[2][j] + numRight[2][j] + 1) <= 5) )
            angleOfArrival[j] = phaseInterferometry(j, locationForAnalysis[j]);
        else
            //printf("SDR1: %d\t", (numLeft[0][j] + numRight[0][j] + 1));
            //printf("SDR2: %d\t", (numLeft[1][j] + numRight[1][j] + 1));
            //printf("SDR3: %d\n", (numLeft[2][j] + numRight[2][j] + 1));
            angleOfArrival[j] = NAN;
    }
    return 0;
}

void setPeakZero(void)
{
    int i, j, k;
    for(i = 0; i < 3; i++){
        for(j = 0; j < NUM_BANDS; j++){
            for(k = 0; k < numFFTs; k++){
                peaks[i][j][k] = 0;
            }
        }
    }
}

void fft_init(void)
{
    int i;

    //Allocate Memory for Data from SDR 1
    fft1n = SAMPLE_TIME;
    fft1in = (fftw_complex *)fftw_malloc(SAMPLE_LENGTH * sizeof(*fft1in));
    fft1out = (fftw_complex *)fftw_malloc(SAMPLE_LENGTH * sizeof(*fft1out));
    for(i = 0; i < SAMPLE_LENGTH; i++){
        fft1in[i] = fft1out[i] = 0;
    }

    //Allocate Memory for Data from SDR 2
    fft2n = SAMPLE_TIME;
    fft2in = (fftw_complex *)fftw_malloc(SAMPLE_LENGTH * sizeof(*fft1in));
    fft2out = (fftw_complex *)fftw_malloc(SAMPLE_LENGTH * sizeof(*fft1out));
    for(i = 0; i < SAMPLE_LENGTH; i++){
        fft2in[i] = fft2out[i] = 0;
    }

    //Allocate Memory for Data from SDR 3
    fft3in = (fftw_complex *)fftw_malloc(SAMPLE_LENGTH * sizeof(*fft1in));
    fft3out = (fftw_complex *)fftw_malloc(SAMPLE_LENGTH * sizeof(*fft1out));
    for(i = 0; i < SAMPLE_LENGTH; i++){
        fft3in[i] = fft3out[i] = 0;
    }

    //Create Plans for All 3 Data Sets
    fft1plan = fftw_plan_many_dft(
                                    1, &fft1n, numFFTs,
                                   fft1in, NULL, 1, distance,
                                   fft1out, NULL, 1, distance,
                                   FFTW_FORWARD, FFTW_ESTIMATE
                                   );

    fft2plan = fftw_plan_many_dft(
                                   1, &fft2n, numFFTs,
                                   fft2in, NULL, 1, distance,
                                   fft2out, NULL, 1, distance,
                                   FFTW_FORWARD, FFTW_ESTIMATE
                                   );
    fft3plan = fftw_plan_many_dft(
                                   1, &fft2n, numFFTs,
                                   fft3in, NULL, 1, distance,
                                   fft3out, NULL, 1, distance,
                                   FFTW_FORWARD, FFTW_ESTIMATE
                                   );
    setPeakZero();
}

int DSP(uint8_t *SDR1_data, uint8_t *SDR2_data, uint8_t *SDR3_data)
{
    int i, j, k;

    uint8_t *buf1 = SDR1_data;
    fftw_complex *floatbuf1 = fft1in;
    j = 0;
    for(i = 0; i < SAMPLE_LENGTH; i++)
        floatbuf1[i] = (buf1[j++] + I*buf1[j++]) - (127.4f+127.4f*I);
    //

    fftw_execute(fft1plan);

    uint8_t *buf2 = SDR2_data;
    fftw_complex *floatbuf2 = fft2in;
    j = 0;
    for(i = 0; i < SAMPLE_LENGTH; i++)
        floatbuf2[i] = (buf2[j++] + I*buf2[j++]) - (127.4f+127.4f*I);
    fftw_execute(fft2plan);

    uint8_t *buf3 = SDR3_data;
    fftw_complex *floatbuf3 = fft3in;
    j = 0;
    for(i = 0; i < SAMPLE_LENGTH; i++)
        floatbuf3[i] = (buf3[j++] + I*buf3[j++]) - (127.4f+127.4f*I);
    fftw_execute(fft3plan);

    //Process the outputs
    //Loop through every FFT taken
    for(i = 0; i < numFFTs; i++){
        //Find the peak of every bin
        for(j = 0; j < distance; j++){
            for(k = 0; k < 3; k++){
                float peak1mag = 0,
                      peak2mag = 0,
                      peak3mag = 0,
                      peak4mag = 0,
                      peak5mag = 0,
                      peak6mag = 0,
                      peak7mag = 0,
                      peak8mag = 0,
                      peak9mag = 0,
                      peak10mag = 0;
                if(k == 0){
                    float complex c = fft1out[i*j + j];
                    float mag = sqrt(crealf(c)*crealf(c) + cimagf(c)*cimagf(c));
                    if((peak1mag < mag) && j < distance/10.0)
                        peaks[k][0][i] = c;
                    if(peak2mag < mag && j < distance/5.0 && j > distance/10.0)
                        peaks[k][1][i] = c;
                    if(peak3mag < mag && j < 3.0*distance/10.0 && j > distance/5.0)
                        peaks[k][2][i] = c;
                    if(peak4mag < mag && j < 4.0*distance/10.0 && j > 3.0*distance/10.0)
                        peaks[k][3][i] = c;
                    if(peak5mag < mag && j < 5.0*distance/10.0 && j > 4.0*distance/10.0)
                        peaks[k][4][i] = c;
                    if(peak6mag < mag && j < 6.0*distance/10.0 && j > 5.0*distance/10.0)
                        peaks[k][5][i] = c;
                    if(peak7mag < mag && j < 7.0*distance/10.0 && j > 6.0*distance/10.0)
                        peaks[k][6][i] = c;
                    if(peak8mag < mag && j < 8.0*distance/10.0 && j > 7.0*distance/10.0)
                        peaks[k][7][i] = c;
                    if(peak9mag < mag && j < 9.0*distance/10.0 && j > 8.0*distance/10.0)
                        peaks[k][8][i] = c;
                    if(peak10mag < mag && j < distance && j > 9.0*distance/10.0)
                        peaks[k][9][i] = c;
                }
                else if(k == 1){
                    float complex c = fft2out[i*j + j];
                    float mag = sqrt(crealf(c)*crealf(c) + cimagf(c)*cimagf(c));
                    if((peak1mag < mag) && j < distance/10.0)
                        peaks[k][0][i] = c;
                    if(peak2mag < mag && j < distance/5.0 && j > distance/10.0)
                        peaks[k][1][i] = c;
                    if(peak3mag < mag && j < 3.0*distance/10.0 && j > distance/5.0)
                        peaks[k][2][i] = c;
                    if(peak4mag < mag && j < 4.0*distance/10.0 && j > 3.0*distance/10.0)
                        peaks[k][3][i] = c;
                    if(peak5mag < mag && j < 5.0*distance/10.0 && j > 4.0*distance/10.0)
                        peaks[k][4][i] = c;
                    if(peak6mag < mag && j < 6.0*distance/10.0 && j > 5.0*distance/10.0)
                        peaks[k][5][i] = c;
                    if(peak7mag < mag && j < 7.0*distance/10.0 && j > 6.0*distance/10.0)
                        peaks[k][6][i] = c;
                    if(peak8mag < mag && j < 8.0*distance/10.0 && j > 7.0*distance/10.0)
                        peaks[k][7][i] = c;
                    if(peak9mag < mag && j < 9.0*distance/10.0 && j > 8.0*distance/10.0)
                        peaks[k][8][i] = c;
                    if(peak10mag < mag && j < distance && j > 9.0*distance/10.0)
                        peaks[k][9][i] = c;
                }
                else if(k == 2){
                    float complex c = fft3out[i*j + j];
                    float mag = sqrt(crealf(c)*crealf(c) + cimagf(c)*cimagf(c));
                    if((peak1mag < mag) && j < distance/10.0)
                        peaks[k][0][i] = c;
                    if(peak2mag < mag && j < distance/5.0 && j > distance/10.0)
                        peaks[k][1][i] = c;
                    if(peak3mag < mag && j < 3.0*distance/10.0 && j > distance/5.0)
                        peaks[k][2][i] = c;
                    if(peak4mag < mag && j < 4.0*distance/10.0 && j > 3.0*distance/10.0)
                        peaks[k][3][i] = c;
                    if(peak5mag < mag && j < 5.0*distance/10.0 && j > 4.0*distance/10.0)
                        peaks[k][4][i] = c;
                    if(peak6mag < mag && j < 6.0*distance/10.0 && j > 5.0*distance/10.0)
                        peaks[k][5][i] = c;
                    if(peak7mag < mag && j < 7.0*distance/10.0 && j > 6.0*distance/10.0)
                        peaks[k][6][i] = c;
                    if(peak8mag < mag && j < 8.0*distance/10.0 && j > 7.0*distance/10.0)
                        peaks[k][7][i] = c;
                    if(peak9mag < mag && j < 9.0*distance/10.0 && j > 8.0*distance/10.0)
                        peaks[k][8][i] = c;
                    if(peak10mag < mag && j < distance && j > 9.0*distance/10.0)
                        peaks[k][9][i] = c;
                }
            }
        }
    }
    findSignal();

    return 0;

}

void findPhaseDifference(uint8_t *SDR1_cal, uint8_t *SDR2_cal, uint8_t *SDR3_cal)
{

    fftw_complex *resultAC;
    fftw_complex *resultBC;

    resultAC = fftw_malloc(2*CALIBRATION_LENGTH * sizeof(*resultAC));
    resultBC = fftw_malloc(2*CALIBRATION_LENGTH * sizeof(*resultBC));
    int i,j;

    uint8_t *buf1 = SDR1_cal;
    fftw_complex *floatbuf1;
    floatbuf1 = fftw_malloc( CALIBRATION_LENGTH * sizeof(*floatbuf1));
    j = 0;
    for(i = 0; i < CALIBRATION_LENGTH; i++)
        floatbuf1[i] = (buf1[j++] + I*buf1[j++]) - (127.4f+127.4f*I);

    uint8_t *buf2 = SDR2_cal;
    fftw_complex *floatbuf2;
    floatbuf2 = fftw_malloc(CALIBRATION_LENGTH * sizeof(*floatbuf2));
    j = 0;
    for(i = 0; i < CALIBRATION_LENGTH; i++){
        floatbuf2[i] = (buf2[j++] + I*buf2[j++]) - (127.4f+127.4f*I);
    }

    uint8_t *buf3 = SDR3_cal;
    fftw_complex *floatbuf3;
    floatbuf3 = fftw_malloc(CALIBRATION_LENGTH * sizeof(*floatbuf3));
    j = 0;
    for(i = 0; i < CALIBRATION_LENGTH; i++)
        floatbuf3[i] = (buf3[j++] + I*buf3[j++]) - (127.4f+127.4f*I);


    xcorr(floatbuf1, floatbuf3, resultAC, CALIBRATION_LENGTH);
    xcorr(floatbuf2, floatbuf3, resultBC, CALIBRATION_LENGTH);

    float maxAC = 0.0, maxBC = 0.0;
    int locAC, locBC;

    for(i = 0; i < (2 * CALIBRATION_LENGTH); i++){
        float magAC = sqrt((crealf(resultAC[i]) * crealf(resultAC[i])) + (cimagf(resultAC[i]) * cimagf(resultAC[i])));
        float magBC = sqrt(crealf(resultBC[i]) * crealf(resultBC[i]) + cimagf(resultBC[i]) * cimagf(resultBC[i]));
        if(magAC > maxAC){
            locAC = i;
            maxAC = magAC;
        }
        if(magBC > maxBC){
            locBC = i;
            maxBC = magBC;
        }
    }

    //printf("AC Loc: %d\n", locAC);
    //printf("BC Loc: %d\n", locBC);
    //printf("Real: %f\n", crealf(resultAC[locAC]));
    //printf("Imag: %f\n", cimagf(resultAC[locAC]));
    phaseCorrectionAC = atan(cimagf(resultAC[locAC])/crealf(resultAC[locAC]))*(180/M_PI);
    //printf("Phase Correction AC before: %f\n", phaseCorrectionAC);
    if(cimagf(resultAC[locAC]) > 0 && crealf(resultAC[locAC]) < 0)
        phaseCorrectionAC = (phaseCorrectionAC * -1.0) + 90.0;
    else if(cimagf(resultAC[locAC]) < 0 && crealf(resultAC[locAC]) < 0)
        phaseCorrectionAC = (phaseCorrectionAC * -1.0) - 90.0;

    phaseCorrectionBC = atan(cimagf(resultBC[locBC])/crealf(resultBC[locBC]))*(180/M_PI);
    //printf("Phase Correction BC before: %f\n", phaseCorrectionBC);
    if(cimagf(resultBC[locBC]) > 0 && crealf(resultBC[locBC]) < 0)
        phaseCorrectionBC = (phaseCorrectionBC * -1.0) + 90.0;
    else if(cimagf(resultBC[locBC]) < 0 && crealf(resultBC[locBC]) < 0)
        phaseCorrectionBC = (phaseCorrectionBC * -1.0) - 90.0;

}
/*
int main(){
    //Open the Data file for testing
    FILE *file;
    int x, i = 0;
    uint8_t       *sdr0Data,
                  *sdr1Data,
                  *sdr2Data;

    sdr0Data = malloc((2*SAMPLE_LENGTH)* sizeof(*sdr0Data));
    sdr1Data = malloc((2*SAMPLE_LENGTH)* sizeof(*sdr0Data));
    sdr2Data = malloc((2*SAMPLE_LENGTH)* sizeof(*sdr0Data));

    fft_init();
    correlateInit();

    file = fopen("sdr0_lastOneV76.dat", "r");

    if (!file){
        printf("There was a problem oping up data set 0\n");
         return 1;
    }
    //printf("The first one is: %hhu\n", sdr0Data[1000]);
    //printf("I am reading the file\n");
    while (fread(sdr0Data, 2*SAMPLE_LENGTH, 1, file));
        //printf("%s",sdr0Data);
        //i++;

    //printf("The first one is: %hhu\n", sdr0Data[1000]);
    fclose(file);

    file = fopen("sdr1_lastOneV76.dat", "r");

    if (!file){
        printf("There was a problem oping up data set 1\n");
        return 1;
    }

    while (fread(sdr1Data, 2*SAMPLE_LENGTH, 1, file));
    //printf("%s",sdr1Data);

    fclose(file);

    file = fopen("sdr2_lastOneV76.dat", "r");

    if (!file){
        printf("There was a problem oping up data set 2\n");
        return 1;
    }

    while (fread(sdr2Data, 2 * SAMPLE_LENGTH, 1, file));
        //printf("%s",sdr2Data);

    fclose(file);
    findPhaseDifference(sdr0Data, sdr1Data, sdr2Data);
    printf("Phase Correction AC: %f\n", phaseCorrectionAC);
    printf("Phase Correction BC: %f\n", phaseCorrectionBC);
    DSP(sdr0Data, sdr1Data, sdr2Data);

    for(i = 0; i < NUM_BANDS; i++)
        printf("Band %d: %f\n", i, angleOfArrival[i]);

    return 0;
}
*/
