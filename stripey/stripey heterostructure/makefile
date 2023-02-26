 CC=g++ -O3
 LD=g++ -O3
 CFLAGS=-O3 -pedantic
 
%.o: %.cpp parameter.h
	 $(CC) -c $(CFLAGS) $<
	 
ALL: conjugategradientperiodicy		bogoliubovdiagonalisation	thermalconductivity

conjugategradientperiodicy: conjugategradientperiodicy.o
		$(LD) $(LDFLAGS) -o conjugategradientperiodicy conjugategradientperiodicy.o -lgsl -lgslcblas -lm
		
bogoliubovdiagonalisation: bogoliubovdiagonalisation.o
		$(LD) $(LDFLAGS) -o bogoliubovdiagonalisation bogoliubovdiagonalisation.o -lgsl -lgslcblas -lm

thermalconductivity: thermalconductivity.o
		$(LD) $(LDFLAGS) -o thermalconductivity thermalconductivity.o -lgsl -lgslcblas -lm
