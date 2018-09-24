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
#define fftLength (SAMPLE_LENGTH/SAMPLE_TIME)

static int nreceivers = 0, corrlen = 0, fft1n = 0, fft2n = 0;

static int sync_debug = 1;

static fftw_complex *fft1in, *fft1out, *fft2in, *fft2out, *fft3in, *fft3out, *fft4in, *fft4out, *fft5in, *fft5out, *fft6in, *fft6out;
static fftw_plan fft1plan, fft2plan, fft3plan, fft4plan, fft5plan, fft6plan;
float  complex peaks[3][NUM_BANDS][numFFTs];
int    peakLocations[3][NUM_BANDS][numFFTs], signalLocation[NUM_BANDS], signalLength = 6*distance;
float  angleOfArrival[NUM_BANDS], phaseCorrectionAC, phaseCorrectionBC, phaseCorrectionAB, phaseAC, phaseBC;

int numberToSubtract = 1;

/* Initialize Correlation */
void correlateInit(){
    phaseCorrectionAC = 0.0;
    phaseCorrectionBC = 0.0;
    phaseCorrectionAB = 0.0;
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

void findPhaseDifferenceSignals(uint8_t *SDR1_cal, uint8_t *SDR2_cal, uint8_t *SDR3_cal, int fftStart){
    fftw_complex *resultAC;
    fftw_complex *resultBC;

    resultAC = fftw_malloc(2*signalLength * sizeof(*resultAC));
    resultBC = fftw_malloc(2*signalLength * sizeof(*resultBC));
    int i,j;

    uint8_t *buf1 = SDR1_cal;
    fftw_complex *floatbuf1;
    floatbuf1 = fftw_malloc( signalLength * sizeof(*floatbuf1));
    j = 0;
    for(i = 0; i < signalLength; i++){
        floatbuf1[i] = (buf1[fftStart + j] + I*buf1[fftStart + j+1]) - (127.4f+127.4f*I);
        j += 2;
    }

    uint8_t *buf2 = SDR2_cal;
    fftw_complex *floatbuf2;
    floatbuf2 = fftw_malloc(signalLength * sizeof(*floatbuf2));
    j = 0;
    for(i = 0; i < signalLength; i++){
        floatbuf2[i] = (buf2[fftStart + j] + I*buf2[fftStart + j+1]) - (127.4f+127.4f*I);
        j += 2;
    }

    uint8_t *buf3 = SDR3_cal;
    fftw_complex *floatbuf3;
    floatbuf3 = fftw_malloc(signalLength * sizeof(*floatbuf3));
    j = 0;
    for(i = 0; i < signalLength; i++){
        floatbuf3[i] = (buf3[fftStart + j] + I*buf3[fftStart + j+1]) - (127.4f+127.4f*I);
        j += 2;
    }


    xcorr(floatbuf1, floatbuf3, resultAC, signalLength);
    xcorr(floatbuf2, floatbuf3, resultBC, signalLength);

    float maxAC = 0.0, maxBC = 0.0;
    int locAC, locBC;

    for(i = 0; i < (2 * signalLength); i++){
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

    phaseAC = atan(cimagf(resultAC[locAC])/crealf(resultAC[locAC]))*(180/M_PI);
    if(cimagf(resultAC[locAC]) > 0 && crealf(resultAC[locAC]) < 0)
        phaseAC = (phaseAC * -1.0) + 90.0;
    else if(cimagf(resultAC[locAC]) < 0 && crealf(resultAC[locAC]) < 0)
        phaseAC = (phaseAC * -1.0) - 90.0;

    phaseBC = atan(cimagf(resultBC[locBC])/crealf(resultBC[locBC]))*(180/M_PI);
    if(cimagf(resultBC[locBC]) > 0 && crealf(resultBC[locBC]) < 0)
        phaseBC = (phaseBC * -1.0) + 90.0;
    else if(cimagf(resultBC[locBC]) < 0 && crealf(resultBC[locBC]) < 0)
        phaseBC = (phaseBC * -1.0) - 90.0;
}

float phaseInterferometry(int band, float *phase){
    int i;
    float top, bottom, phaseAC, phaseBC;
    //top = 1.1547 * (phase[0] - phase[2] + phaseCorrectionAC);
    //bottom = phase[1] - phase[2] + phaseCorrectionBC;
    printf("Phase CorrectionAC: %f\t Phase Correction BC: %f\n", phaseCorrectionAC, phaseCorrectionBC);
    phaseAC = phase[0] - phase[2] + phaseCorrectionAC;
    phaseBC = phase[1] - phase[2] + phaseCorrectionBC;
printf("Phase AC: %f\t Phase BC: %f\n", phaseAC, phaseBC);
    //if((top > 0 && bottom > 0) || (top < 0 && bottom < 0))
    //    offset = 30.0;

    if(phaseAC >360.0)
        phaseAC =(phaseAC - 360.0);
    else if(phaseAC < -360.0)
        phaseAC += 360.0;
    if(phaseBC > 360.0)
        phaseBC -= 360.0;
    else if(phaseBC < -360.0)
        phaseBC += 360.0;

    top = (1.40963 * phaseAC) - (0.70482 * phaseBC);
    bottom = 1.22078 * phaseBC;

    printf("Phase 0: %f\t Phase 1: %f\t Phase 2: %f\n", phase[0], phase[1], phase[2]);
    printf("Before Top is: %f\n",top);
    printf("Before Bottom is: %f\n",bottom);

    /*if(top > 180.0)
        top -= 360.0;
    else if(top < -180.0)
        top += 360.0;
    if(bottom > 180.0)
        bottom -= 360.0;
    else if(bottom < -180.0)
        bottom += 360.0;*/

    /*if(top > 180){
        int topMult = top/180;
        top -= topMult*180;
    }
    else if(top < 180){
        int topMult = top/(-180);
        top += topMult*180;
    }

    if(bottom > 180){
        int botMult = bottom/180;
        bottom -= botMult*180;
    }
    else if(bottom < 180){
        int botMult = bottom/(-180);
        bottom += botMult*180;
    }*/


    //top = top * 1.1547;

    /*if(top > 90.0)
        top -= 180.0;
    else if(top < -90.0)
        top += 180.0;

    if(bottom > 90.0)
        bottom -= 180.0;
    else if(bottom < -90.0)
        bottom += 180.0;*/

    printf("Top is: %f\n",top);
    printf("Bottom is: %f\n",bottom);

    /*if(bottom < 0 && top < 0)
        offset = -180.0;
    else if(bottom > 0 && top < 0)
        offset = 180.0;*/

    return atanf(top/bottom)*(180/M_PI);

}



/*float phaseInterferometry(int band, int *location){
    int i;
    float phase[3], top, bottom;
    for(i = 0; i < 3; i++){
        phase[i] = (atanf(cimagf(peaks[i][band][location[i]])/crealf(peaks[i][band][location[i]]))) * (180.0/M_PI);
        if(cimagf(peaks[i][band][location[i]]) > 0 && crealf(peaks[i][band][location[i]]) < 0)
            phase[i] = (phase[i] * -1.0) + 90.0;
        else if(cimagf(peaks[i][band][location[i]]) < 0 && crealf(peaks[i][band][location[i]]) < 0)
            phase[i] = (phase[i] * -1.0) - 90.0;
    }
    top = 1.1547 * (phase[0] - phase[2] + phaseCorrectionAC);
    bottom = phase[1] - phase[2] + phaseCorrectionBC;

    printf("Top is: %f\n",top);
    printf("Bottom is: %f\n",bottom);

    if(top > 180.0)
        top -= 360.0;
    else if(top < -180.0)
        top += 360.0;
    if(bottom > 180.0)
        bottom -= 360.0;
    else if(bottom < -180.0)
        bottom += 360.0;

    printf("Top is: %f\n",top);
    printf("Bottom is: %f\n",bottom);

    return atanf((top/bottom) - 0.57735)*(180/M_PI);

}*/

void signalFFT(uint8_t *SDR1_data, uint8_t *SDR2_data, uint8_t *SDR3_data, int fftStart){

    int i, j;
    //int signalLength = 8 * distance;

    //Does not do ffts if no signal
    if(fftStart == -1)
        return;
    printf("fftStart: %d\n", fftStart);
    uint8_t *buf1 = SDR1_data;
    fftw_complex *floatbuf1 = fft4in;
    j = 0;
    for(i = 0; i < signalLength; i++){
        floatbuf1[i] = (buf1[2*fftStart + j] + I*buf1[2*fftStart + j+1]) - (127.4f+127.4f*I);
        j += 2;
    }
    fftw_execute(fft4plan);

    uint8_t *buf2 = SDR2_data;
    fftw_complex *floatbuf2 = fft5in;
    j = 0;
    for(i = 0; i < signalLength; i++){
        floatbuf2[i] = (buf2[2*fftStart + j] + I*buf2[2*fftStart + j+1]) - (127.4f+127.4f*I);
        j += 2;
    }
    fftw_execute(fft5plan);

    uint8_t *buf3 = SDR3_data;
    fftw_complex *floatbuf3 = fft6in;
    j = 0;
    for(i = 0; i < signalLength; i++){
        floatbuf3[i] = (buf3[2*fftStart + j] + I*buf3[2*fftStart + j+1]) - (127.4f+127.4f*I);
        j += 2;
    }
    fftw_execute(fft6plan);

}

void findSignalPeaks(int band){
    float complex c;
    double mag;
    int lower = (COLLAR_OFFSET-COLLAR_TOLERANCE)/(SAMPLE_FREQUENCY/(signalLength/2));
    int upper = (COLLAR_OFFSET+COLLAR_TOLERANCE)/(SAMPLE_FREQUENCY/(signalLength/2));
    int start = signalLength*(band)/24 + lower, stop = signalLength*(band)/24 + upper;
    //int start = signalLength*(band-1)/24, stop = signalLength;
    int i, j, k;
    float peak = 0;
    float complex finalSignal[3];
    int finalSignalLocation[3];

    //printf("Signal Length %d\t Start: %d\t Stop: %d\t Center: %d\n", signalLength, start, stop, signalLength*band/24);
    for(k = 0; k < 3; k++){
        peak = 0;
        for(j = start; j < stop; j++){
            if(k == 0)
                c = fft4out[j];
            else if(k == 1)
                c = fft5out[j];
            else if(k == 2)
                c = fft6out[j];
            mag = sqrt(crealf(c)*crealf(c) + cimagf(c)*cimagf(c));
            if(peak < mag){
                peak = mag;
                finalSignal[k] = c;
                finalSignalLocation[k] = j;
            }
        }

    }

    /*if(finalSignalLocation[0] <= finalSignalLocation[1] && finalSignalLocation[0] <= finalSignalLocation[2])
        start = finalSignalLocation[0];
    else if(finalSignalLocation[1] <= finalSignalLocation[0] && finalSignalLocation[1] <= finalSignalLocation[2])
        start = finalSignalLocation[1];
    else if(finalSignalLocation[2] <= finalSignalLocation[1] && finalSignalLocation[2] <= finalSignalLocation[0])
        start = finalSignalLocation[2];*/

    /*if(finalSignalLocation[0] >= finalSignalLocation[1] && finalSignalLocation[0] >= finalSignalLocation[2])
        stop = finalSignalLocation[0];
    else if(finalSignalLocation[1] >= finalSignalLocation[0] && finalSignalLocation[1] >= finalSignalLocation[2])
        stop = finalSignalLocation[1];
    else if(finalSignalLocation[2] >= finalSignalLocation[1] && finalSignalLocation[2] >= finalSignalLocation[0])
        stop = finalSignalLocation[2];*/

    peak = 0;
    float mag1, mag2, mag3;
    float complex c1, c2, c3;
    int finalSignalLocationTemp;
    for(j = start; j < stop+1; j++){
            c1 = fft4out[j];
            c2 = fft5out[j];
            c3 = fft6out[j];
        mag1 = sqrt(crealf(c1)*crealf(c1) + cimagf(c1)*cimagf(c1));
        mag2 = sqrt(crealf(c2)*crealf(c2) + cimagf(c2)*cimagf(c2));
        mag3 = sqrt(crealf(c3)*crealf(c3) + cimagf(c3)*cimagf(c3));
        mag = (mag1 + mag2 + mag3)/3;
        if(peak < mag){
            peak = mag;
            finalSignalLocationTemp = j;
            //printf("%d Mag: %f\n", j, mag);
        }
        /*if(j == 19131)
            printf("19131 Mag: %f\n", mag);
        else if(j == 19147 || j == 19148 || j == 19184 || j == 19183)
            printf("%d Mag: %f\n", j, mag);*/
    }
    printf("Location Temp %d\n", finalSignalLocationTemp);

    if(finalSignalLocation[0] == finalSignalLocation[1] && finalSignalLocation[0] != finalSignalLocation[2])
        finalSignalLocation[2] =finalSignalLocation[0];
    else if(finalSignalLocation[0] == finalSignalLocation[2] && finalSignalLocation[0] != finalSignalLocation[1])
        finalSignalLocation[1] = finalSignalLocation[0];
    else if(finalSignalLocation[2] == finalSignalLocation[1] && finalSignalLocation[2] != finalSignalLocation[0])
        finalSignalLocation[0] = finalSignalLocation[2];

    finalSignal[0] = fft4out[finalSignalLocationTemp];
    finalSignal[1] = fft5out[finalSignalLocationTemp];
    finalSignal[2] = fft6out[finalSignalLocationTemp];
    /*finalSignal[0] = fft4out[4139];
    finalSignal[1] = fft5out[4139];
    finalSignal[2] = fft6out[4139];*/
    for(k = 0; k < 3; k++)
        printf("Location %d: %d\t Peak: %f\n", k, finalSignalLocation[k], sqrt(crealf(finalSignal[k])*crealf(finalSignal[k]) + cimagf(finalSignal[k])*cimagf(finalSignal[k])));
    //if(finalSignalLocation[0] == finalSignalLocation[1] && finalSignalLocation[0] == finalSignalLocation[2]){
        float phase[3];
        for(i = 0; i < 3; i++){
            phase[i] = (atanf(cimagf(finalSignal[i])/crealf(finalSignal[i]))) * (180.0/M_PI);
            if(cimagf(finalSignal[i]) > 0 && crealf(finalSignal[i]) < 0)
                phase[i] = 180.0 + phase[i];
            else if(cimagf(finalSignal[i]) < 0 && crealf(finalSignal[i]) < 0)
                phase[i] = (180.0 - phase[i]) * -1.0;
            //printf("Phase %d: %f\n", i, phase[i]);
        }
        angleOfArrival[band] = phaseInterferometry(band, phase);
    //}

}

float findSignal(){
    //Initializ Variables
    int i, j, k, maxLocation[3][NUM_BANDS] = {{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}};
    float max[3][NUM_BANDS] = {{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}};
    float averageNoise[3][NUM_BANDS] = {{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}};

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


    int numLeft[3][NUM_BANDS], numRight[3][NUM_BANDS];
    int badLeftFlag = 0, badRightFlag = 0, stopLeft = 0, stopRight = 0;

    //Now we verify that the signal exists
    for(i = 0; i < 3; i++){
        for(j = 0; j < NUM_BANDS; j++){
            stopLeft = 0;
            stopRight = 0;
            numLeft[i][j] = 0;
            numRight[i][j] = 0;
            badRightFlag = 0;
            badLeftFlag = 0;
            float threshold = max[i][j] - ( (2*max[i][j]-averageNoise[i][j]) /3 );
            for(k = 1; k < 4; k++){
                float peakmagMin = sqrt(crealf(peaks[i][j][maxLocation[i][j] - k])*crealf(peaks[i][j][maxLocation[i][j] - k]) + cimagf(peaks[i][j][maxLocation[i][j] - k])*cimagf(peaks[i][j][maxLocation[i][j] - k]));
                float peakmagPlus = sqrt(crealf(peaks[i][j][maxLocation[i][j] + k])*crealf(peaks[i][j][maxLocation[i][j] + k]) + cimagf(peaks[i][j][maxLocation[i][j] + k])*cimagf(peaks[i][j][maxLocation[i][j] + k]));

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

    int location[3][NUM_BANDS], finalNumLeft[3][NUM_BANDS], finalNumRight[3][NUM_BANDS];
    k = 0;
    int numberRight = 0, numberLeft = 0;
    for(i = 0; i < 3; i++){
        for(j = 0; j < NUM_BANDS; j++){
            int curLocationLeft = peakLocations[i][j][maxLocation[i][j]];
            int curLocationRight = peakLocations[i][j][maxLocation[i][j]];
            while(curLocationLeft == peakLocations[i][j][maxLocation[i][j]]){
                curLocationLeft = peakLocations[i][j][maxLocation[i][j] - k];
                numberLeft = k;
                k++;
            }
            k = 0;
            while(curLocationRight == peakLocations[i][j][maxLocation[i][j]]){
                curLocationRight = peakLocations[i][j][maxLocation[i][j] + k];
                numberRight = k;
                k++;
            }
            finalNumRight[i][j] = numberRight;
            finalNumLeft[i][j] = numberLeft;
            location[i][j] = numberRight + numberLeft;
           // if(location[i][j] <= 5 && location[i][j] >= 3)
           //     printf("Band %d: SDR%d Left: %d Right: %d Location: %d Max: %d\n", j, i, numberLeft, numberRight, peakLocations[i][j][maxLocation[i][j]], maxLocation[i][j] );
        }
    }
    float complex phaseComplex[3][NUM_BANDS];
    float phase[3][NUM_BANDS];
    for(i = 0; i < NUM_BANDS; i++){
        phaseComplex[0][i] = fft1out[maxLocation[0][i]*distance + peakLocations[0][i][maxLocation[0][i]]];
        phaseComplex[1][i] = fft2out[maxLocation[0][i]*distance + peakLocations[0][i][maxLocation[0][i]]];
        phaseComplex[2][i] = fft3out[maxLocation[0][i]*distance + peakLocations[0][i][maxLocation[0][i]]];

        phase[0][i] = (atanf(cimagf(phaseComplex[0][i])/crealf(phaseComplex[0][i]))) * (180.0/M_PI);
        if(cimagf(phaseComplex[0][i]) > 0 && crealf(phaseComplex[0][i]) < 0)
            phase[0][i] = (phase[0][i] * -1.0) + 90;
        else if(cimagf(phaseComplex[0][i]) < 0 && crealf(phaseComplex[0][i]) < 0)
            phase[0][i] = (phase[0][i] * -1.0) - 90.0;

        phase[1][i] = (atanf(cimagf(phaseComplex[1][i])/crealf(phaseComplex[1][i]))) * (180.0/M_PI);
        if(cimagf(phaseComplex[1][i]) > 0 && crealf(phaseComplex[1][i]) < 0)
            phase[1][i] = (phase[1][i] * -1.0) + 90.0;
        else if(cimagf(phaseComplex[1][i]) < 0 && crealf(phaseComplex[1][i]) < 0)
            phase[1][i] = (phase[1][i] * -1.0) - 90.0;

        phase[2][i] = (atanf(cimagf(phaseComplex[2][i])/crealf(phaseComplex[2][i]))) * (180.0/M_PI);
        if(cimagf(phaseComplex[2][i]) > 0 && crealf(phaseComplex[2][i]) < 0)
            phase[2][i] = (phase[2][i] * -1.0) + 90.0;
        else if(cimagf(phaseComplex[2][i]) < 0 && crealf(phaseComplex[2][i]) < 0)
            phase[2][i] = (phase[2][i] * -1.0) - 90.0;
    }

    //Pick the points to use
    int locationForAnalysis[NUM_BANDS];


    for(j = 0; j < NUM_BANDS; j++){
        if(((maxLocation[0][j] - numLeft[0][j]) >= (maxLocation[1][j] - numLeft[1][j])) && ((maxLocation[0][j] - numLeft[0][j]) >= (maxLocation[2][j] - numLeft[2][j])))
            locationForAnalysis[j] = maxLocation[0][j] - finalNumLeft[0][j]-numberToSubtract;
        else if(((maxLocation[1][j] - numLeft[1][j]) >= (maxLocation[0][j] - numLeft[0][j])) && ((maxLocation[1][j] - numLeft[1][j]) >= (maxLocation[2][j] - numLeft[2][j])))
            locationForAnalysis[j] = (maxLocation[1][j] - finalNumLeft[1][j])-numberToSubtract;
        else if(((maxLocation[2][j] - numLeft[2][j]) >= (maxLocation[0][j] - numLeft[0][j])) && ((maxLocation[2][j] - numLeft[2][j]) >= (maxLocation[1][j] - numLeft[1][j])))
            locationForAnalysis[j] = (maxLocation[2][j] - finalNumLeft[2][j])-numberToSubtract;

        if( ((numLeft[0][j] + numRight[0][j] + 1) >= 3) && ((numLeft[0][j] + numRight[0][j] + 1) <= 5) && ((numLeft[1][j] + numRight[1][j] + 1) >= 3) && ((numLeft[1][j] + numRight[1][j] + 1) <= 5) &&((numLeft[2][j] + numRight[2][j] + 1) >= 3) && ((numLeft[2][j] + numRight[2][j] + 1) <= 5) && (location[0][j] >= 4 && location[0][j] <= 5) && (location[1][j] >= 4 && location[1][j] <= 5) && (location[2][j] >=4 && location[0][j] <=5)){
            int location[3] = {(maxLocation[0][j] - finalNumLeft[0][j] + 1),(maxLocation[1][j] - finalNumLeft[1][j] + 1),(maxLocation[2][j] - finalNumLeft[2][j] + 1)};
            signalLocation[j] = distance * locationForAnalysis[j];
            printf("SDR1 Max: %f\t SDR2 Max: %f\t SDR3 Max: %f\n", max[0][j], max[1][j], max[2][j]);
            printf("Location for Analysis: %d\n", locationForAnalysis[j]);
        }
        else{
            angleOfArrival[j] = NAN;
            signalLocation[j] = -1;
        }
    }
    return 0;
}

void setPeakZero(){
    int i, j, k;
    for(i = 0; i < 3; i++){
        for(j = 0; j < NUM_BANDS; j++){
            for(k = 0; k < numFFTs; k++){
                peaks[i][j][k] = 0;
                peakLocations[i][j][k] = 0;
            }
        }
    }
}

void fft_init(){
    int i;
    //int signalLength = 8*distance;

    //Allocate Memory for Data from SDR 1
    fft1n = distance;
    fft1in = (fftw_complex *)fftw_malloc(SAMPLE_LENGTH * sizeof(*fft1in));
    fft1out = (fftw_complex *)fftw_malloc(SAMPLE_LENGTH * sizeof(*fft1out));
    for(i = 0; i < SAMPLE_LENGTH; i++){
        fft1in[i] = fft1out[i] = 0;
    }

    //Allocate Memory for Data from SDR 2
    fft2n = distance;
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

    fft4in = (fftw_complex *)fftw_malloc(signalLength * sizeof(*fft1in));
    fft4out = (fftw_complex *)fftw_malloc(signalLength * sizeof(*fft1out));
    for(i = 0; i < signalLength; i++){
        fft4in[i] = fft4out[i] = 0;
    }

    fft5in = (fftw_complex *)fftw_malloc(signalLength * sizeof(*fft1in));
    fft5out = (fftw_complex *)fftw_malloc(signalLength * sizeof(*fft1out));
    for(i = 0; i < signalLength; i++){
        fft5in[i] = fft5out[i] = 0;
    }

    fft6in = (fftw_complex *)fftw_malloc(signalLength * sizeof(*fft1in));
    fft6out = (fftw_complex *)fftw_malloc(signalLength * sizeof(*fft1out));
    for(i = 0; i < signalLength; i++){
        fft6in[i] = fft6out[i] = 0;
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
    fft4plan = fftw_plan_dft_1d(signalLength, fft4in, fft4out, FFTW_FORWARD, FFTW_ESTIMATE);
    fft5plan = fftw_plan_dft_1d(signalLength, fft5in, fft5out, FFTW_FORWARD, FFTW_ESTIMATE);
    fft6plan = fftw_plan_dft_1d(signalLength, fft6in, fft6out, FFTW_FORWARD, FFTW_ESTIMATE);
    setPeakZero();
}

void findFFTPeaks(int sdrNumber, int band, int fftNum){
    float complex c;
    float mag, peakMag = 0.0;
    int i, j, k = sdrNumber;

    int lower = (COLLAR_OFFSET-COLLAR_TOLERANCE)/(SAMPLE_FREQUENCY/(distance/2));
    int upper = (COLLAR_OFFSET+COLLAR_TOLERANCE)/(SAMPLE_FREQUENCY/(distance/2));
    int start = distance*(band)/24 + lower, stop = distance*(band)/24 + upper;
    //int start = distance*(band)/24, stop = distance*(band+1)/24;

    for(j = start; j < stop; j++){
        if(sdrNumber == 0)
            c = fft1out[fftNum*distance + j];
        else if(sdrNumber == 1)
            c = fft2out[fftNum*distance + j];
        else if(sdrNumber == 2)
            c = fft3out[fftNum*distance + j];
        mag = sqrt(crealf(c)*crealf(c) + cimagf(c)*cimagf(c));
        if((peakMag < mag)){
            peakMag = mag;
            peaks[k][band][fftNum] = c;
            peakLocations[k][band][fftNum] = j;
        }
    }
}

int DSP(uint8_t *SDR1_data, uint8_t *SDR2_data, uint8_t *SDR3_data){
    int i, j, k;

    uint8_t *buf1 = SDR1_data;
    fftw_complex *floatbuf1 = fft1in;
    j = 0;
    for(i = 0; i < SAMPLE_LENGTH; i++){
        floatbuf1[i] = (buf1[j] + I*buf1[j+1]) - (127.4f+127.4f*I);
        j += 2;
    }
    fftw_execute(fft1plan);

    uint8_t *buf2 = SDR2_data;
    fftw_complex *floatbuf2 = fft2in;
    j = 0;
    for(i = 0; i < SAMPLE_LENGTH; i++){
        floatbuf2[i] = (buf2[j] + I*buf2[j+1]) - (127.4f+127.4f*I);
        j += 2;
    }
    fftw_execute(fft2plan);

    uint8_t *buf3 = SDR3_data;
    fftw_complex *floatbuf3 = fft3in;
    j = 0;
    for(i = 0; i < SAMPLE_LENGTH; i++){
        floatbuf3[i] = (buf3[j] + I*buf3[j+1]) - (127.4f+127.4f*I);
        j += 2;
    }
    fftw_execute(fft3plan);

    //Process the outputs
    //Loop through every FFT taken

    for(i = 0; i < 3; i++){
        for(j = 0; j < numFFTs; j++){
            for(k = 0; k < NUM_BANDS; k++){
                //Find the peak of every bin
                findFFTPeaks(i, k, j);
            }
        }
    }
    findSignal();

    for(i = 0; i < NUM_BANDS; i++){
        signalFFT(SDR1_data, SDR2_data, SDR3_data, signalLocation[i]);
        if(signalLocation[i] != -1){
            findSignalPeaks(i);
            //findPhaseDifferenceSignals(SDR1_data, SDR2_data, SDR3_data, signalLocation[i]);
            //printf("Phase Difference AC: %f\t Phase Difference BC: %f\n", phaseAC, phaseBC);
        }
    }

    fftw_free(fft1in);
    fftw_free(fft2in);
    fftw_free(fft3in);
    fftw_free(fft1out);
    fftw_free(fft2out);
    fftw_free(fft3out);

    fftw_free(fft4in);
    fftw_free(fft5in);
    fftw_free(fft6in);
    fftw_free(fft4out);
    fftw_free(fft5out);
    fftw_free(fft6out);

    for(i = 0; i < NUM_BANDS; i++)
        printf("Band %d: %f\n", i, angleOfArrival[i]);


    return 0;

}

void findPhaseDifference(uint8_t *SDR1_cal, uint8_t *SDR2_cal, uint8_t *SDR3_cal){

    fftw_complex *resultAC;
    fftw_complex *resultBC;
    fftw_complex *resultAB;

    resultAC = fftw_malloc(2*CALIBRATION_LENGTH * sizeof(*resultAC));
    resultBC = fftw_malloc(2*CALIBRATION_LENGTH * sizeof(*resultBC));
    resultAB = fftw_malloc(2*CALIBRATION_LENGTH * sizeof(*resultAB));
    int i,j;

    uint8_t *buf1 = SDR1_cal;
    fftw_complex *floatbuf1;
    floatbuf1 = fftw_malloc( CALIBRATION_LENGTH * sizeof(*floatbuf1));
    j = 0;
    for(i = 0; i < CALIBRATION_LENGTH; i++){
        floatbuf1[i] = (buf1[j] + I*buf1[j+1]) - (127.4f+127.4f*I);
        j += 2;
    }

    uint8_t *buf2 = SDR2_cal;
    fftw_complex *floatbuf2;
    floatbuf2 = fftw_malloc(CALIBRATION_LENGTH * sizeof(*floatbuf2));
    j = 0;
    for(i = 0; i < CALIBRATION_LENGTH; i++){
        floatbuf2[i] = (buf2[j] + I*buf2[j+1]) - (127.4f+127.4f*I);
        j += 2;
    }

    uint8_t *buf3 = SDR3_cal;
    fftw_complex *floatbuf3;
    floatbuf3 = fftw_malloc(CALIBRATION_LENGTH * sizeof(*floatbuf3));
    j = 0;
    for(i = 0; i < CALIBRATION_LENGTH; i++){
        floatbuf3[i] = (buf3[j] + I*buf3[j+1]) - (127.4f+127.4f*I);
        j += 2;
    }


    xcorr(floatbuf1, floatbuf3, resultAC, CALIBRATION_LENGTH);
    xcorr(floatbuf2, floatbuf3, resultBC, CALIBRATION_LENGTH);
    xcorr(floatbuf1, floatbuf2, resultAB, CALIBRATION_LENGTH);

    float maxAC = 0.0, maxBC = 0.0, maxAB = 0.0;
    int locAC, locBC, locAB;

    for(i = 0; i < (2 * CALIBRATION_LENGTH); i++){
        float magAC = sqrt((crealf(resultAC[i]) * crealf(resultAC[i])) + (cimagf(resultAC[i]) * cimagf(resultAC[i])));
        float magBC = sqrt(crealf(resultBC[i]) * crealf(resultBC[i]) + cimagf(resultBC[i]) * cimagf(resultBC[i]));
        float magAB = sqrt((crealf(resultAB[i]) * crealf(resultAB[i])) + (cimagf(resultAB[i]) * cimagf(resultAB[i])));
        if(magAC > maxAC){
            locAC = i;
            maxAC = magAC;
        }
        if(magBC > maxBC){
            locBC = i;
            maxBC = magBC;
        }
        if(magAB > maxAB){
            locAB = i;
            maxAB = magAB;
        }
    }
   
    phaseCorrectionAC = atan(cimagf(resultAC[locAC])/crealf(resultAC[locAC]))*(180/M_PI);
    printf("Phase Correction AC Before: %f\n", phaseCorrectionAC);
    if(cimagf(resultAC[locAC]) > 0 && crealf(resultAC[locAC]) < 0)
        phaseCorrectionAC = 180.0 + phaseCorrectionAC;
    else if(cimagf(resultAC[locAC]) < 0 && crealf(resultAC[locAC]) < 0)
        phaseCorrectionAC = (180.0 - phaseCorrectionAC) * -1.0;

    phaseCorrectionBC = atan(cimagf(resultBC[locBC])/crealf(resultBC[locBC]))*(180/M_PI);
    if(cimagf(resultBC[locBC]) > 0 && crealf(resultBC[locBC]) < 0)
        phaseCorrectionBC = (180.0 + phaseCorrectionBC);
    else if(cimagf(resultBC[locBC]) < 0 && crealf(resultBC[locBC]) < 0)
        phaseCorrectionBC = (180.0 - phaseCorrectionBC) * -1.0;

    phaseCorrectionAB = atan(cimagf(resultAB[locAB])/crealf(resultAB[locAB]))*(180/M_PI);
    if(cimagf(resultAB[locAB]) > 0 && crealf(resultAB[locAB]) < 0)
        phaseCorrectionAB = (180.0 + phaseCorrectionAB);
    else if(cimagf(resultAB[locAB]) < 0 && crealf(resultAB[locAB]) < 0)
        phaseCorrectionAB = (180.0 - phaseCorrectionAB) * -1.0;
	
    fftw_free(resultAC);
    fftw_free(resultAB);
    fftw_free(resultBC);
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
    while (fread(sdr0Data, 2*SAMPLE_LENGTH, 1, file));

    fclose(file);

    file = fopen("sdr1_lastOneV76.dat", "r");

    if (!file){
        printf("There was a problem oping up data set 1\n");
        return 1;
    }

    while (fread(sdr1Data, 2*SAMPLE_LENGTH, 1, file));

    fclose(file);

    file = fopen("sdr2_lastOneV76.dat", "r");

    if (!file){
        printf("There was a problem oping up data set 2\n");
        return 1;
    }

    while (fread(sdr2Data, 2 * SAMPLE_LENGTH, 1, file));

    fclose(file);
    findPhaseDifference(sdr0Data, sdr1Data, sdr2Data);
    printf("Phase Correction AC: %f\n", phaseCorrectionAC);
    printf("Phase Correction BC: %f\n", phaseCorrectionBC);
    printf("Phase Correction AB: %f\n", phaseCorrectionAB);
    DSP(sdr0Data, sdr1Data, sdr2Data);


    return 0;
}
*/
