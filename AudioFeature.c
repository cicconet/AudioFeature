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

#include "AudioFeature.h"

// ----------------------------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------------------------

double * AFWaveRead(const char * fileName, int * numFrames)
{
	int blockSize = 512;
	
	SNDFILE	* infile = NULL ;
	SF_INFO	sfinfo ;
	
	if ((infile = sf_open(fileName, SFM_READ, &sfinfo)) == NULL)	{
		printf("Not able to open input file %s.\n", fileName);
		puts(sf_strerror(NULL));
		return NULL ;
	} ;
	
	printf("frames: %ld\n", (long int)sfinfo.frames);
	printf("samplerate: %d\n", sfinfo.samplerate);
	printf("channels: %d\n", sfinfo.channels);
	printf("format: %d\n", sfinfo.format);
	printf("sections: %d\n", sfinfo.sections);
	printf("seekable: %d\n", sfinfo.seekable);
	
	// read data
	double * data = (double *)calloc(sfinfo.frames, sizeof(double));
	
	double buf[sfinfo.channels*blockSize];
	int k, m, readcount;
	double nChannels = (double)sfinfo.channels;
	double value;
	long int frame = 0;
	while ((readcount = sf_readf_double(infile, buf, blockSize)) > 0) {
		for (k = 0 ; k < readcount ; k++) {
			value = 0;
			for (m = 0 ; m < sfinfo.channels ; m++) {
				value += (buf[k*sfinfo.channels+m]/nChannels);
			}
			data[frame] = value;
			frame += 1;
		}
	}
	
	sf_close (infile) ;
	
	*numFrames = sfinfo.frames;
	return data;
}

AFBuffer * AFAllocScratchBuffer(int windowSize)
{
	AFBuffer * buffer = (AFBuffer *)malloc(sizeof(AFBuffer));
	buffer->windowSize = windowSize;
	buffer->halfWindowSize = windowSize/2;
	buffer->complexVector = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*windowSize);
	buffer->hannWindow = (double *)malloc(windowSize*sizeof(double));
	// buffer->realVectorFullSize = (double *)malloc(windowSize*sizeof(double));
	buffer->realVectorHalfSize = (double *)malloc(buffer->halfWindowSize*sizeof(double));
	
	double center = (double)buffer->halfWindowSize;
	double width = (double)windowSize;
	for (int i = 0; i < windowSize; i++) {
		double point = (double)i;
		buffer->hannWindow[i] = evaluateHannWindow(center, width, point);
	}
	
	return buffer;
}

void AFReleaseScratchBuffer(AFBuffer * scratchBuffer)
{
	free(scratchBuffer->hannWindow);
	// free(scratchBuffer->realVectorFullSize);
	free(scratchBuffer->realVectorHalfSize);
	fftw_free(scratchBuffer->complexVector);
	free(scratchBuffer);
}

void AFFFTMagnitudes(AFBuffer * scratchBuffer, double * magnitudes, double * audioWindow)
{
	in_place_forward_fft_real(scratchBuffer, audioWindow);
	magnitudes_of_complex_vector(magnitudes, scratchBuffer->complexVector, scratchBuffer->halfWindowSize);
}

void AFPowerSpectrum(AFBuffer * scratchBuffer, double * powerSpectrum, double * audioWindow)
{
	in_place_forward_fft_real(scratchBuffer, audioWindow);
	sq_magnitudes_of_complex_vector(powerSpectrum, scratchBuffer->complexVector, scratchBuffer->halfWindowSize);
}

double AFLoudness(AFBuffer * scratchBuffer, double * audioWindow)
{
	in_place_forward_fft_real(scratchBuffer, audioWindow);
	sq_magnitudes_of_complex_vector(scratchBuffer->realVectorHalfSize, scratchBuffer->complexVector, scratchBuffer->halfWindowSize);
	return sum(scratchBuffer->realVectorHalfSize, scratchBuffer->halfWindowSize);
}

double AFSpectralCentroid(AFBuffer * scratchBuffer, double * audioWindow)
{
	in_place_forward_fft_real(scratchBuffer, audioWindow);
	magnitudes_of_complex_vector(scratchBuffer->realVectorHalfSize, scratchBuffer->complexVector, scratchBuffer->halfWindowSize);
	double s = sum(scratchBuffer->realVectorHalfSize, scratchBuffer->halfWindowSize);
	double specCentBin = 0.0;
	if (s > 0.0) {
		for (int i = 0; i < scratchBuffer->halfWindowSize; i++) {
			specCentBin += (double)(i+1)*scratchBuffer->realVectorHalfSize[i]/s;
		}
	}
	scratchBuffer->sumFFTMag = s;
	return specCentBin;
}

double AFSpectralSpread(AFBuffer * scratchBuffer, double * audioWindow)
{
	double specCentBin = AFSpectralCentroid(scratchBuffer, audioWindow);
	double s = scratchBuffer->sumFFTMag;
	double specSprdBin = 0.0;
	double difference;
	if (s > 0.0) {
		for (int i = 0; i < scratchBuffer->halfWindowSize; i++) {
			difference = (double)i+1.0-specCentBin;
			difference *= difference;
			specSprdBin += difference*scratchBuffer->realVectorHalfSize[i]/s;
		}
	}
	return specSprdBin;
}

double AFSpectralSkewness(AFBuffer * scratchBuffer, double * audioWindow)
{
	double specCentBin = AFSpectralCentroid(scratchBuffer, audioWindow);
	double s = scratchBuffer->sumFFTMag;
	double specSprdBin = 0.0;
	double difference;
	if (s > 0.0) {
		for (int i = 0; i < scratchBuffer->halfWindowSize; i++) {
			difference = (double)i+1.0-specCentBin;
			difference *= difference;
			specSprdBin += difference*scratchBuffer->realVectorHalfSize[i]/s;
		}
	}
	double sigma = sqrt(specSprdBin);
	double specSknsBin = 0.0;
	if (s > 0.0) {
		for (int i = 0; i < scratchBuffer->halfWindowSize; i++) {
			difference = (double)i+1.0-specCentBin;
			difference = difference*difference*difference;
			specSknsBin += difference*scratchBuffer->realVectorHalfSize[i]/s;
		}
		specSknsBin = specSknsBin/(sigma*sigma*sigma);
	}
	return specSknsBin;
}

double AFSpectralRollOff(AFBuffer * scratchBuffer, double * audioWindow)
{
	double loudness = AFLoudness(scratchBuffer, audioWindow);
	double spectralRollOff = 0.0;
	if (loudness > 0.0) {
		double rollOff = 0.95*loudness;
		double accumulator = 0.0;
		int bin = 0;
		while (accumulator < rollOff) {
			accumulator += scratchBuffer->realVectorHalfSize[bin];
			bin += 1;
		}
		spectralRollOff = (double)bin;
	}
	return spectralRollOff;
}

double AFSpectralFlatness(AFBuffer * scratchBuffer, double * audioWindow)
{
	double s, p, specFlatBin;
	double exp = 1.0/(double)scratchBuffer->halfWindowSize;
	in_place_forward_fft_real(scratchBuffer, audioWindow);
	magnitudes_of_complex_vector(scratchBuffer->realVectorHalfSize, scratchBuffer->complexVector, scratchBuffer->halfWindowSize);
	s = sum(scratchBuffer->realVectorHalfSize, scratchBuffer->halfWindowSize);
	specFlatBin = 0.0;
	if (s > 0.0) {
		p = 1.0;
		for (int i = 0; i < scratchBuffer->halfWindowSize; i++) {
			p *= pow(scratchBuffer->realVectorHalfSize[i], exp);
		}
		specFlatBin = scratchBuffer->halfWindowSize*p/s;
	}
	return specFlatBin;
}

double AFSpectralCrest(AFBuffer * scratchBuffer, double * audioWindow)
{
	double s, m, specCrestBin;
	in_place_forward_fft_real(scratchBuffer, audioWindow);
	magnitudes_of_complex_vector(scratchBuffer->realVectorHalfSize, scratchBuffer->complexVector, scratchBuffer->halfWindowSize);
	s = sum(scratchBuffer->realVectorHalfSize, scratchBuffer->halfWindowSize);
	m = max(scratchBuffer->realVectorHalfSize, scratchBuffer->halfWindowSize);
	specCrestBin = 0.0;
	if (s > 0.0) {
		specCrestBin = scratchBuffer->halfWindowSize*m/s;
	}
	return specCrestBin;
}

void AFPreChroma(AFBuffer * scratchBuffer, double * preChroma, double * audioWindow)
{
	// WARNING: halfWindowSize possible values: 512, 1024, 2048, ...

	AFPowerSpectrum(scratchBuffer, scratchBuffer->realVectorHalfSize, audioWindow);

	double frequencyLag = 22050.0/(double)scratchBuffer->halfWindowSize;
	int halfNumIndBellowWindow = floor(32.703194*scratchBuffer->halfWindowSize/22050.0)+1;
	int numIndBellowWindow = 2*halfNumIndBellowWindow;
	double hannWindowSize = (numIndBellowWindow-1)*frequencyLag;
	
	// compute the 84 piano frequencies probability vector (a.k.a. pre-chroma vector)
	double frequency, accumulator, hannValue, currentFrequency;
	double sevenOctavesDensity[84]; // 7*12 = 84
	double m = 0.0;
	int currentIndex, toTheLeftOfFreqIndex;
	for (int i = -45; i < 39; i++) { // 84 piano frequency values
		frequency = 440.0*pow(2.0, i/12.0); // equivalent midi note: 69+i
		toTheLeftOfFreqIndex = floor(frequency/frequencyLag);
		accumulator = 0.0;
		for (int j = 0; j < numIndBellowWindow; j++) {
			currentIndex = toTheLeftOfFreqIndex-halfNumIndBellowWindow+1+j;
			currentFrequency = currentIndex*frequencyLag;
			hannValue = evaluateHannWindow(frequency, hannWindowSize, currentFrequency);
			accumulator = accumulator+hannValue*scratchBuffer->realVectorHalfSize[currentIndex];
		}
		sevenOctavesDensity[i+45] = accumulator;
		if (accumulator > m) m = accumulator;
	}
	
	// normalize pre chroma vector
	if (m > 0.0) {
		for (int i = 0; i < 84; i++) {
			preChroma[i] = sevenOctavesDensity[i]/m;
		}
	}
	else {
		for (int i = 0; i < 84; i++) {
			preChroma[i] = 0.0;
		}
	}
}

void AFChroma(AFBuffer * scratchBuffer, double * chroma, double * audioWindow)
{
	AFPowerSpectrum(scratchBuffer, scratchBuffer->realVectorHalfSize, audioWindow);
	
	int i = -45;
	double minFrequency = 440.0*powf(2.0, (double)(i-0.5)/12.0);
	i = 75;
	double maxFrequency = 440.0*powf(2.0, ((double)i-0.5)/12.0);
	
	double frequency, nearestExactFrequency, distance;
	double fIndex;
	int index;
	int binsPerFrequencyInterval[120]; for (int i = 0; i < 120; i++) { binsPerFrequencyInterval[i] = 0; }
	double preChromaVector[120]; for (int i = 0; i < 120; i++) { preChromaVector[i] = 0.0; }
	for (int i = 1; i < scratchBuffer->halfWindowSize; i++) {
		frequency = 22050.0*(double)i/(double)scratchBuffer->halfWindowSize;
		if (frequency >= minFrequency && frequency < maxFrequency) {
			fIndex = log2(pow(frequency/440.0, 12.0))+0.5;
			index = (int)floor(fIndex);
			nearestExactFrequency = 440.0*pow(2.0, (double)index/12.0);
			distance = abs(frequency-nearestExactFrequency);
			binsPerFrequencyInterval[index+45] += 1;
			preChromaVector[index+45] += chromaWeightForDistance(distance)*scratchBuffer->realVectorHalfSize[i];
		}
	}
	for (int i = 0; i < 120; i++) {
		if (binsPerFrequencyInterval[i] > 0) {
			preChromaVector[i] /= (double)binsPerFrequencyInterval[i];
		} else {
			preChromaVector[i] = 0.0;
		}
	}
	float m = 0.0;
	for (int key = 0; key < 12; key++) {
		chroma[key] = 0.0;
		for (int scale = 0; scale < 10; scale++) {
			chroma[key] += preChromaVector[scale*12+key];
		}
		if (chroma[key] > m) m = chroma[key];
	}
	if (m > 0) {
		for (int key = 0; key < 12; key++) {
			chroma[key] /= m;
		}
	}
}

// ----------------------------------------------------------------------------------------------------
// Auxiliary
// ----------------------------------------------------------------------------------------------------

double chromaWeightForDistance(double distance)
{
	return exp(-distance*distance/1000.0);
}

double evaluateHannWindow(double center, double width, double point)
{
	return (1.0-cos(6.28318531*(point-center+width/2.0)/width))/width;
}

double sum(double * vector, int length)
{
	double s = 0.0;
	for (int i = 0; i < length; i++) {
		s += vector[i];
	}
	return s;
}

double max(double * vector, int length)
{
	double m = -INFINITY;
	for (int i = 0; i < length; i++) {
		if (vector[i] > m) {
			m = vector[i];
		}
	}
	return m;
}

void sq_magnitudes_of_complex_vector(double * sq_magnitudes, fftw_complex * complexVector, int length)
{
	for (int i = 0; i < length; i++) {
		double x = complexVector[i][0];
		double y = complexVector[i][1];
		sq_magnitudes[i] = x*x+y*y;
	}	
}

void magnitudes_of_complex_vector(double * magnitudes, fftw_complex * complexVector, int length)
{
	for (int i = 0; i < length; i++) {
		double x = complexVector[i][0];
		double y = complexVector[i][1];
		magnitudes[i] = sqrt(x*x+y*y);
	}
}

void in_place_forward_fft_real(AFBuffer * scratchBuffer, double * audioWindow)
{
	copy_real(scratchBuffer->complexVector, audioWindow, scratchBuffer->windowSize, 1, scratchBuffer->hannWindow);
	clean_im(scratchBuffer->complexVector, scratchBuffer->windowSize);
	in_place_forward_fft_complex(scratchBuffer->complexVector, scratchBuffer->windowSize);
}

void in_place_forward_fft_complex(fftw_complex * complexVector, int length)
{
	fftw_plan p = fftw_plan_dft_1d(length, complexVector, complexVector, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
}

void clean_im(fftw_complex * complex_array, int length)
{
	for (int i = 0; i < length; i++) {
		complex_array[i][1] = 0.0;
	}
}

void copy_real(fftw_complex * complex_array, double * input, int length, int windowing, double * window)
{
	if (windowing) {
		for (int i = 0; i < length; i++) {
			complex_array[i][0] = input[i]*window[i];
		}
	} else {
		for (int i = 0; i < length; i++) {
			complex_array[i][0] = input[i];
		}
	}
}

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
