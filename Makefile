CXX=g++ 
CXXFLAGS=-g -O2 -std=c++0x -ftree-vectorize -msse2 -ftree-vectorizer-verbose=1 -ffast-math 
INCLUDE=-I/home/bkandel/workspace/gen3-fftw-fit-2 
FITSLIBS=-lCCfits -lcfitsio
GSLLIBS=-lgsl -lgslcblas
FFTWLIBS=-lfftw3

gen3-fftw-fit-2: gen3-fftw-fit-2.cpp paramlist.cc paramlist.h  
	${CXX} ${CXXFLAGS} gen3-fftw-fit-2.cpp paramlist.cc -o gen3-fftw-fit-2.o ${INCLUDE} ${FITSLIBS} ${FFTWLIBS}
gen3-fftw-fit-2-init: gen3-fftw-fit-2-init.cpp paramlist.cc paramlist.h
	${CXX} ${CXXFLAGS} gen3-fftw-fit-2-init.cpp paramlist.cc -o gen3-fftw-ffit-2-init.o ${INCLUDE} ${FITSLIBS} ${FFTWLIBS}

