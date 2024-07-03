## To compile, select the compile options below and enter 'make' into the terminal
## the command 'make clean' removes all compiled files and executable files for recompilation

# for gfortran Compiler
#======================
F90          = gfortran
F90LINKER    = gfortran
DEFS      =
#FFLAGS   	= -Og -pipe -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all -Wconversion -fbounds-check -ffpe-trap=invalid -llapack
FFLAGS   = -O3 -pipe -llapack
INCLUDES  =
LFLAGS    = $(FFLAGS)

# for nvfortran Compiler
#======================
# F90          = nvfortran
# F90LINKER    = nvfortran
# DEFS      =
# FFLAGS   = -fast
# INCLUDES  =
# LFLAGS    = $(FFLAGS)


# for ifort Compiler
#======================
# F90          = ifort
# F90LINKER    = ifort
# DEFS      =
# FFLAGS   = -O3
# INCLUDES  =
# LFLAGS    = $(FFLAGS)


OBJECTS = \
	  WENO4_mod.o \
	  call_twostr.o \
	  k_Rosseland_mod.o \
	  IC_mod.o \
	  dry_conv_adj_mod.o \
		MLT_mod.o \
		sw_direct_mod.o \
		sw_adding_mod.o \
		sw_SDA_mod.o \
		sw_Toon_mod.o \
		sw_SH2_mod.o \
		sw_SH4_mod.o \
		sw_disort_ts_mod.o \
		lw_AA_E_mod.o \
		lw_AA_L_mod.o \
		lw_sc_linear_mod.o \
		lw_VIM_mod.o \
		lw_Toon_mod.o \
		lw_disort_ts_mod.o \
		lw_Feautrier_mod.o \
		lw_DFE_mod.o \
    FMS_RC.o

# executable statement
EXECS  = ../FMS_RC_sg

.SUFFIXES: .o .f90 .F .f

default: FMS_RC_sg

FMS_RC_sg:  $(OBJECTS)
	 $(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)

clean:
	rm -f *.o *.mod *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.F.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

