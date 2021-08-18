
# Specify Location of HDF5 binaries (with include/ and lib/ underneath)
HDF5=/opt/hdf5

# Compiler
# Use the HDF Fortran compiler wrapper....
# FC = gfortran
FC = $(HDF5)/bin/h5fc

FFLAGS = -O2 -fbacktrace -Wall -fbounds-check

.SUFFIXES :
.SUFFIXES : .f90 .o .exe

EXEC = makewb.exe

default: $(EXEC)

$(EXEC) : Mod_hdftools.o Mod_geom.o Mod_state.o main.o writehdf.o
	$(FC) $(FFLAGS) $^ -o $@

Mod_hdftools.o : ../src/Mod_hdftools.f90
	$(FC) $(FFLAGS) -c $^

.f90.o :
	$(FC) $(FFLAGS) -c $*.f90

clean:
	rm -f *.mod *.o $(EXEC)

main.o : Mod_geom.o Mod_state.o
writehdf.o : Mod_hdftools.o Mod_geom.o Mod_state.o

