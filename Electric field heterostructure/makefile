 CC=g++ -O3
 LD=g++ -O3
 CFLAGS=-O3 -pedantic
 
%.o: %.cpp parameter.h
	 $(CC) -c $(CFLAGS) $<
	 
ALL: conjugategradient	conjugategradientperiodicy	angleplot	gradient	bogoliubovdiagonalisation	thermalconductivity	spintexture

conjugategradient: conjugategradient.o
		$(LD) $(LDFLAGS) -o conjugategradient conjugategradient.o -lgsl -lgslcblas -lm

conjugategradientperiodicy: conjugategradientperiodicy.o
		$(LD) $(LDFLAGS) -o conjugategradientperiodicy conjugategradientperiodicy.o -lgsl -lgslcblas -lm

gradient: gradient.o
		$(LD) $(LDFLAGS) -o gradient gradient.o -lgsl -lgslcblas -lm

angleplot: angleplot.o
		$(LD) $(LDFLAGS) -o angleplot angleplot.o -lgsl -lgslcblas -lm
		
bogoliubovdiagonalisation: bogoliubovdiagonalisation.o
		$(LD) $(LDFLAGS) -o bogoliubovdiagonalisation bogoliubovdiagonalisation.o -lgsl -lgslcblas -lm

thermalconductivity: thermalconductivity.o
		$(LD) $(LDFLAGS) -o thermalconductivity thermalconductivity.o -lgsl -lgslcblas -lm

spintexture: spintexture.o
		$(LD) $(LDFLAGS) -o spintexture spintexture.o -lgsl -lgslcblas -lm
		
