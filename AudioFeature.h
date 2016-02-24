/*
 *  AudioFeature Library
 *	Compute some low-level audio descriptors.
 *  
 *  Marcelo Cicconet, Visgraf / IMPA
 *	www.impa.br/~cicconet
 *
 *	This library is distributed under the GNU General Public Licence Version 3.
 *	See notice at the end of this file.
 *
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sndfile.h>
#include <math.h>
#include <fftw3.h>

struct AFBuffer {
	int windowSize;
	int halfWindowSize;
	double sumFFTMag;
	fftw_complex * complexVector;
	double * hannWindow;
	// double * realVectorFullSize;
	double * realVectorHalfSize;
};
typedef struct AFBuffer AFBuffer;

// ----------------------------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------------------------
double * AFWaveRead(const char * fileName, int * numFrames);
// read audio file with extension .wav (returns the sum of the channels)

AFBuffer * AFAllocScratchBuffer(int windowSize);
// buffer needed for all subsequent computations

void AFReleaseScratchBuffer(AFBuffer * scratchBuffer);
// frees scratchBuffer allocated space

void AFFFTMagnitudes(AFBuffer * scratchBuffer, double * magnitudes, double * audioWindow);
// returns the magnitudes of the DFT of audioWindow; length of magnitudes is length(audioWindow)/2

void AFPowerSpectrum(AFBuffer * scratchBuffer, double * powerSpectrum, double * audioWindow);
// returns the squared magnitudes of the DFT of audioWindow; length of powerSpectrum is length(audioWindow)/2

double AFLoudness(AFBuffer * scratchBuffer, double * audioWindow);
// sum of the squared magnitudes of the spectrum of frequencies

double AFSpectralCentroid(AFBuffer * scratchBuffer, double * audioWindow);
// barycenter of the magnitudes of the spectrum of frequencies

double AFSpectralSpread(AFBuffer * scratchBuffer, double * audioWindow);
// variance of the magnitudes of the spectrum of frequencies

double AFSpectralSkewness(AFBuffer * scratchBuffer, double * audioWindow);
// 3rd order moment of the magnitudes of the spectrum of frequencies
// skewness = 0: symmetric distribution
// skewness < 0: more energy on the right
// skewness > 0: more energy on the left

double AFSpectralRollOff(AFBuffer * scratchBuffer, double * audioWindow);
// frequency so that 95% of the signal energy is contained bellow it
// it is correlated somehow to the noise cutting frequency

double AFSpectralFlatness(AFBuffer * scratchBuffer, double * audioWindow);
// ratio of the geometric mean to the arithmetic mean of the magnitudes of the spectrum of frequencies
// its a measure of the noisiness (flatness) of the spectrum of the spectrum magnitudes

double AFSpectralCrest(AFBuffer * scratchBuffer, double * audioWindow);
// ratio of the maximum value to the arithmetic mean of the magnitudes of the spectrum of frequencies

void AFPreChroma(AFBuffer * scratchBuffer, double * preChroma, double * audioWindow); // 84-dimensional
// 84 piano frequency energies (from midi note 24 to midi note 107)
// entries are in the range [0,1]

void AFChroma(AFBuffer * scratchBuffer, double * chroma, double * audioWindow); // 12-dimensional
// frequency energies for each of the 12 keys (from C to B)
// entries are in the range [0,1]

// ----------------------------------------------------------------------------------------------------
// Auxiliary
// ----------------------------------------------------------------------------------------------------
double chromaWeightForDistance(double distance);
double evaluateHannWindow(double center, double width, double point);
double sum(double * vector, int length);
double max(double * vector, int length);
void sq_magnitudes_of_complex_vector(double * sq_magnitudes, fftw_complex * complexVector, int length);
void magnitudes_of_complex_vector(double * magnitudes, fftw_complex * complexVector, int length);
void in_place_forward_fft_real(AFBuffer * scratchBuffer, double * audioWindow);
void in_place_forward_fft_complex(fftw_complex * complexVector, int length);
void clean_im(fftw_complex * complex_array, int length);
void copy_real(fftw_complex * complex_array, double * input, int length, int windowing, double * window);

/*
 *	Copyright (C) 2011 Marcelo Cicconet
 *
 *	This program is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */