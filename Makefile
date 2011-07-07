CXX=g++ 
CXXFLAGS=-g -O2 -std=c++0x -ftree-vectorize -msse2 -ftree-vectorizer-verbose=1 -ffast-math 
INCLUDE=-I/home/bkandel/include 
FITSLIBS=-lCCfits -lcfitsio
GSLLIBS=-lgsl -lgslcblas
FFTWLIBS=-lfftw3

gen3-fit-fftw: gen3-fit-fftw.cc paramlist.cc paramlist.h 
	${CXX} ${CXXFLAGS} gen3-fit-fftw.cc paramlist.cc -o gen3-fit-fftw ${INCLUDE} ${FITSLIBS} ${FFTWLIBS}

gim-reorient: gim-reorient.cc paramlist.cc paramlist.h
	${CXX} ${CXXFLAGS} gim-reorient.cc paramlist.cc -o gim-reorient ${INCLUDE} ${FITSLIBS} ${GSLLIBS} ${FFTWLIBS} 
 
