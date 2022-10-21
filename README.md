# KQED

This library contains the position-space QED kernel for the Hadronic Light-by-Light contribution to the muon g-2

## Building the code

The library is built with gnu automake, so the typical behaviour is:

    ./configure CC=gcc CFLAGS="-O3 -Wall -fopenmp -mavx -mfma" --prefix={path to install dir}

The code has been checked (compiles without warning under -Wall -Wpedantic) with icc-2018, clang-6, clang-14, gcc-7.4, gcc-8.1, and gcc-11.2. Currently the clang-compiled and the intel-compiled codes run faster.

If using icc I recommend the flag -xHOST instead of -mavx -mfma. The compiler will probably tell you as much. In general I recommend using the icc compiler where available. I strongly recommend building with AVX/FMA if possible as it is significantly faster, likewise with OpenMP.

The code compiles without warning under clang's static analyzer (scan-build), and is memory-leak free according to valgrind.

For the initial build, please run:

       aclocal && autoreconf -ivf

To generate the Makefile.in files.

Also you can enable the loop unrolling using DUFF's device

     --enable-DUFF

although this may be slower than not using is (it seems like with icc this is the case, gcc perhaps not)

And then

    make all install

You will notice that there is no longer an input file as this is all included in the form-factor file which will be copied into your install dir.

The "-mavx" should be supported by your compiler and is turned on to access the crc32c vector instruction for the checksums and the AVX routines in simd.c, most newer processors support the FMA routines as well. All the files being read in are now check-summed for verification.

## Running the code

The file amu_for_lattice.c gives an example of how to use the library and some regression/timing tests. Please after compiling run:

    ./REGress.sh

to make sure your build is sane. For instance when compiling with the flag "--fast-math" in gcc we have seen some very minor round-off differences from our regression files.
