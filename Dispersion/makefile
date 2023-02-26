CC=g++
LD=g++
CFLAGS=-O3 -pedantic

%.o: %.cpp
	$(CC) -c $(CFLAGS) $<

ALL: dispersion plotdispersion

dispersion: dispersion.o
	    $(LD) $(LDFLAGS) -o dispersion dispersion.o -lgsl -lgslcblas -lm

plotdispersion: plotdispersion.o
	    $(LD) $(LDFLAGS) -o plotdispersion plotdispersion.o -lgsl -lgslcblas -lm
