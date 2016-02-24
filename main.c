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

#include <stdio.h>
#include "AudioFeature.h"

int	main (int argc, char * argv [])
{	
	char * fileName = "/Users/Cicconet/Desktop/Song.wav";
	
	double * audioData;
	int numFrames;
	audioData = AFWaveRead(fileName, &numFrames);
	
	int windowSize = 1024;
	int halfWindowSize = windowSize/2;
	int hopSize = halfWindowSize;
	int numWindows = floorf((float)(numFrames-windowSize)/(float)hopSize);
	
	AFBuffer * afBuffer = AFAllocScratchBuffer(windowSize);

	FILE * f = fopen("/Users/Cicconet/Desktop/Features.txt", "w");
	for (int window = 0; window < numWindows; window++) {
		fprintf(f, "%lf\n", AFLoudness(afBuffer, &audioData[window*hopSize]));
	}
	fclose(f);
	
	AFReleaseScratchBuffer(afBuffer);
	
	free(audioData);
	return 0 ;
}

/*
 Copyright (C) 2011 Marcelo Cicconet
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
