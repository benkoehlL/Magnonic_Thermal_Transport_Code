 CC=g++ -O3
 LD=g++ -O3
 CFLAGS=-O3 -pedantic
 
%.o: %.cpp parameter.h
	 $(CC) -c $(CFLAGS) $<
	 
ALL: bogoliubovdiagonalisation

bogoliubovdiagonalisation: bogoliubovdiagonalisation.o
		$(LD) $(LDFLAGS) -o bogoliubovdiagonalisation bogoliubovdiagonalisation.o -lgsl -lgslcblas -lm
