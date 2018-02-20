# Makefile for h2ep. T.Pope, Feb 2018
FC          = gfortran 
FC_ASIS     = $(FC) 
FLAGS_DEBUG = 
FLAGS_OPT   = -O3 
FLAGS_OMP   = -fopenmp 
FLAGS_FORM  = -ffixed-line-length-none 
LAPACK      = /usr/lib/liblapack.so.3
BLAS        = /usr/lib/libblas.so.3 -fexternal-blas
LIBRARY     = $(BLAS) $(LAPACK) 

FLAGS       =  $(FLAGS_DEBUG) $(FLAGS_OPT) $(FLAGS_OMP) $(FLAGS_FORM) 

.f.o:
	$(FC) -c $(FLAGS) $(FLAGS_DEBUG) $<

OBJS = modules.o getgreens.o cutuphamiltonian.o desingularize.o getscatteringmatrix.o initialize_lead.o inverteffham.o make_k.o moses.o readhamiltonian.o subspace_iron.o h2ep.o

ai:	$(OBJS)
	$(FC) -o /home/tom/Desktop/h2ep-code-dev $(OBJS) $(FLAGS) $(LIBRARY)

clean:
	rm *.o; rm *.mod
