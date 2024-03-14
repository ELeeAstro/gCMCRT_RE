## To compile, select the compile options below and enter 'make' into the terminal
## the command 'make clean' removes all compiled files and executable files for recompilation

# for nvfortran Compiler
#======================
F90          = nvfortran
F90LINKER    = nvfortran
DEFS      =
#FFLAGS = -O0 -g -traceback -cuda -cudalib=curand
FFLAGS   = -fast -O3 -cuda -cudalib=curand
INCLUDES  =
LFLAGS    = $(FFLAGS)


OBJECTS = \
mc_precision_mod.o \
mc_data_mod.o \
mc_aux_mod.o \
mc_moments_mod.o \
mc_tauint_mod.o \
mc_inc_mod.o \
mc_emit_mod.o \
mc_scat_mod.o \
grey.o \
semi_grey.o \
gCMCRT_RE.o 

# executable statement
EXECS  = ../gCMCRT_RE

.SUFFIXES: .o .f90 .F

default: gCMCRT_RE

gCMCRT_RE:  $(OBJECTS)
	 $(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)

clean:
	rm -f *.o *.mod *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.F.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
