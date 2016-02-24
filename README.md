# audiofeature
A simple C library for computing low-level audio descriptors.

For more details, see http://w3.impa.br/~cicconet/audiofeature/index.html

Dependencies: FFTW library for DFT computation, and Libsndfile for reading the bytes of a WAVE or AIF file.

The implementation assumes that the audio is sampled at 44100 frames per second.

I'm not providing documentation, since I think AudioFeature is small enough for the code to be self explained.
