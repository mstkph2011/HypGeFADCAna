#vpath %.h ../
#vpath %.cxx ../
#vpath %.h $(ROOTSYS)/include

INCLFLAGS  = $(shell root-config --cflags) \
								-I$(ROOTSYS)/include -I../
LDLIBS     = $(shell root-config --glibs) \
                -L/$(ROOTSYS)/lib -lMathCore -lSpectrum 

#LDLIBS     = $(shell root-config --glibs)
# -lMinuit

#MAINFILE = RealSpecAna.cxx
#MAINFILE = main_ExtractRisetimeCorrParametersAndRisetimePlots.cxx
MAINFILE = main_MogonIIJülichAna.cxx

SOURCES = $(wildcard *.cxx)
HEADERS = $(wildcard *.h)
OBJECTS = THypGePeakFitFunction.o THypGeSpectrumAnalyser.o RealSpecAna
# OBJECTS = $(patsubst %.cxx,%.o,$(wildcard *.cxx))

CC = g++
CFLAGS = -O2 -DNDEBUG -Wall


all:                    RealSpecAna
RealSpecAna: $(MAINFILE) THypGeSpectrumAnalyser.o THypGePeakFitFunction.o
	$(CC) $(CFLAGS) $(INCLFLAGS) $(LDLIBS) $(MAINFILE) THypGeSpectrumAnalyser.o THypGePeakFitFunction.o $(LDLIBS) -o RealSpecAna

THypGeSpectrumAnalyser.o: THypGeSpectrumAnalyser.cxx THypGeSpectrumAnalyser.h
	$(CC) $(CFLAGS) $(INCLFLAGS) -c THypGeSpectrumAnalyser.cxx

THypGePeakFitFunction.o: THypGePeakFitFunction.cxx THypGePeakFitFunction.h
	$(CC) $(CFLAGS) $(INCLFLAGS) -c THypGePeakFitFunction.cxx
.PHONY :                clean

clean :
			rm $(OBJECTS)

