RDIST=/software/R-2.9.0/

INCLUDES= -I. -I$(RDIST)/include -I/usr/lib64/R/include -I/usr/share/R/include/
LIBS = -lg2c

PROF= #-pg
DEB=  -g
OPT = -O 

CPPFLAGS= -pthread -fPIC -Wall $(OPT) $(PROF) $(DEB) $(INCLUDES)

fast_multinom_v2.so: fast_multinom_weight_v2.cpp
	R CMD SHLIB fast_multinom_weight_v2.cpp -o fast_multinom_weight_v2.so


read_minimac.so: read_minimac.cpp
	R CMD SHLIB read_minimac.cpp -o read_minimac.so



hmm.so: hmm.cpp
	R CMD SHLIB hmm.cpp -o hmm.so

parse_vcf.so: parse_vcf.cpp
	R CMD SHLIB parse_vcf.cpp -o parse_vcf.so

map_variant.so: map_variants_to_freq.cpp
	R CMD SHLIB map_variants_to_freq.cpp -o map_variants_to_freq.so


routines.so: consecutive.cpp
	R CMD SHLIB consecutive.cpp  write_beagle.cpp -o routines.so

consecutive.o : consecutive.cpp 

write_beagle.o: write_beagle.cpp

clean:
	rm -f routines.so *.o

listing:
	nm -g routines.so
