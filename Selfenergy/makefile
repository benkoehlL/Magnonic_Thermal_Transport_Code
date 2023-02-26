CC=g++ -O3
LD=g++ -O3

CFLAGS=-O3 -pedantic

%.o: %.cpp parameter.h
	$(CC) -c $(CFLAGS) $<

ALL: lookupvertex readvertex uvtest lookupselfenergy plotselfenergy dos interpolation plotinterpolateddecayrate thermalconductivity densityplot reducedthermalconductivity templookupselfenergy tempplotselfenergy tempinterpolation tempdensityplot compactlookupselfenergy plotcompactselfenergy compactinterpolation bfieldthermalconductivity powerlawplot templookupselfenergyconstB tempinterpolationconstB thermalconductivityconstB dispersion plotdispersion

lookupvertex: lookupvertex.o
	    $(LD) $(LDFLAGS) -o lookupvertex lookupvertex.o -lgsl -lgslcblas -lm

readvertex: readvertex.o
	    $(LD) $(LDFLAGS) -o readvertex readvertex.o -lgsl -lgslcblas -lm

uvtest: uvtest.o
	    $(LD) $(LDFLAGS) -o uvtest uvtest.o -lgsl -lgslcblas -lm

lookupselfenergy: lookupselfenergy.o
	    $(LD) $(LDFLAGS) -o lookupselfenergy lookupselfenergy.o -lgsl -lgslcblas -lm
	    
plotselfenergy: plotselfenergy.o
	    $(LD) $(LDFLAGS) -o plotselfenergy plotselfenergy.o -lgsl -lgslcblas -lm
	    
plotinterpolateddecayrate: plotinterpolateddecayrate.o
	    $(LD) $(LDFLAGS) -o plotinterpolateddecayrate plotinterpolateddecayrate.o -lgsl -lgslcblas -lm
	    
dos: dos.o
	    $(LD) $(LDFLAGS) -o dos dos.o -lgsl -lgslcblas -lm
	    
interpolation: interpolation.o
	    $(LD) $(LDFLAGS) -o interpolation interpolation.o -lgsl -lgslcblas -lm
	    
thermalconductivity: thermalconductivity.o
	    $(LD) $(LDFLAGS) -o thermalconductivity thermalconductivity.o -lgsl -lgslcblas -lm
	    
reducedthermalconductivity: reducedthermalconductivity.o
	    $(LD) $(LDFLAGS) -o reducedthermalconductivity reducedthermalconductivity.o -lgsl -lgslcblas -lm

densityplot: densityplot.o
	    $(LD) $(LDFLAGS) -o densityplot densityplot.o -lgsl -lgslcblas -lm
	    
templookupselfenergy: templookupselfenergy.o
	    $(LD) $(LDFLAGS) -o templookupselfenergy templookupselfenergy.o -lgsl -lgslcblas -lm
	    
tempplotselfenergy: tempplotselfenergy.o
	    $(LD) $(LDFLAGS) -o tempplotselfenergy tempplotselfenergy.o -lgsl -lgslcblas -lm	  
	 
tempinterpolation: tempinterpolation.o
	    $(LD) $(LDFLAGS) -o tempinterpolation tempinterpolation.o -lgsl -lgslcblas -lm	 
	    
tempdensityplot: tempdensityplot.o
	    $(LD) $(LDFLAGS) -o tempdensityplot tempdensityplot.o -lgsl -lgslcblas -lm	  
	    
compactlookupselfenergy: compactlookupselfenergy.o
	    $(LD) $(LDFLAGS) -o compactlookupselfenergy compactlookupselfenergy.o -lgsl -lgslcblas -lm	

plotcompactselfenergy: plotcompactselfenergy.o
	    $(LD) $(LDFLAGS) -o plotcompactselfenergy plotcompactselfenergy.o -lgsl -lgslcblas -lm	
	    
compactinterpolation: compactinterpolation.o
	    $(LD) $(LDFLAGS) -o compactinterpolation compactinterpolation.o -lgsl -lgslcblas -lm
	    
bfieldthermalconductivity: bfieldthermalconductivity.o
	    $(LD) $(LDFLAGS) -o bfieldthermalconductivity bfieldthermalconductivity.o -lgsl -lgslcblas -lm

powerlawplot: powerlawplot.o
	    $(LD) $(LDFLAGS) -o powerlawplot powerlawplot.o -lgsl -lgslcblas -lm

templookupselfenergyconstB: templookupselfenergyconstB.o
	    $(LD) $(LDFLAGS) -o templookupselfenergyconstB templookupselfenergyconstB.o -lgsl -lgslcblas -lm

tempinterpolationconstB: tempinterpolationconstB.o
	    $(LD) $(LDFLAGS) -o tempinterpolationconstB tempinterpolationconstB.o -lgsl -lgslcblas -lm

thermalconductivityconstB: thermalconductivityconstB.o
	    $(LD) $(LDFLAGS) -o thermalconductivityconstB thermalconductivityconstB.o -lgsl -lgslcblas -lm

dispersion: dispersion.o
	    $(LD) $(LDFLAGS) -o dispersion dispersion.o -lgsl -lgslcblas -lm

plotdispersion: plotdispersion.o
	    $(LD) $(LDFLAGS) -o plotdispersion plotdispersion.o -lgsl -lgslcblas -lm
    
run:
	./lookupvertex && ./lookupselfenergy && ./interpolation && ./plotselfenergy && ./plotinterpolateddecayrate && ./thermalconductivity && ./reducedthermalconductivity && ./densityplot && ./templookupselfenergy && ./tempplotselfenergy

runthermal:
	./templookupselfenergyconstB && ./tempinterpolationconstB && ./thermalconductivityconstB

clean:
	rm lookupvertex.o readvertex.o uvtest.o lookupselfenergy.o plotselfenergy.o dos.o interpolation.o plotinterpolateddecayrate.o thermalconductivity.o densityplot.o reducedthermalconductivity.o templookupselfenergy.o tempplotselfenergy.o tempinterpolation.o tempdensityplot.o compactlookupselfenergy.o plotcompactselfenergy.o compactinterpolation.o bfieldthermalconductivity.o powerlawplot.o templookupselfenergyconstB.o tempinterpolationconstB.o thermalconductivityconstB.o dispersion.o plotdispersion.o
