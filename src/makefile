# Makefile for program ace
#
# Once compiled, it is convenient to move the executable to a default OS X path, e.g.:
# > sudo su
# > ...enter password
# > cp ace /usr/local/bin
#============================================================

 FC      = gfortran
 FFLAGS  = -O2 -fno-second-underscore -Wall -fbackslash

 BINDIR  = .
 EXEC    = ace

 OBJECTS = ace_data.o \
           ace.o \
           compute_chemical_equilibrium.o \
           minimize_gibbs_energy.o \
           numrecipes.o \
           mnewt_aux.o \
           read_input.o \
           write_results.o

 COMMON  = ace.common \
           dimension.common

all: $(EXEC)
	echo "             Done!"

$(EXEC) : $(OBJECTS) $(COMMON)
	$(FC) $(FFLAGS) $(OBJECTS) -o $(BINDIR)/$(EXEC)

clean:
	rm -f -r *.o *.mod $(EXEC)
