## To compile, select the compile options below and enter 'make' into the terminal
## the command 'make clean' removes all compiled files and executable files for recompilation

# for gfortran Compiler
#======================
F90          = gfortran
F90LINKER    = gfortran
DEFS      =
#FFLAGS   	= -Og -pipe -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all -Wconversion -fbounds-check
FFLAGS   = -O3 -pipe
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
	  toms715.o \
	  linspace.o \
	  interp_mod.o \
          phys.o \
          read_opacity_tables.o \
          read_planck.o \
          cloud_mod.o \
	  partitioner.o \
	  tau_struc.o \
	  k_Rosseland_mod.o \
	  IC_mod.o \
	  dry_conv_adj_mod.o \
	  moist_adjustment.o \
	  call_twostr.o \
	  ts_isothermal_mod.o \
	  ts_isothermal_2_mod.o \
	  ts_Toon_mod.o \
    ts_Toon_scatter_mod.o \
	  ts_Heng_mod.o \
	  ts_short_char_mod_linear.o \
		ts_short_char_mod_Bezier.o \
		ts_short_char_mod_Bezier_bg.o \
    ts_Mendonca_mod.o \
    ts_Lewis_scatter_mod.o \
	  ts_disort_scatter_mod.o \
	  ts_disort_scatter_bg_mod.o \
    FMS_RC.o

#		ts_Heng_ITS_mod.o \


# executable statement
EXECS  = FMS_RC

.SUFFIXES: .o .f90 .F

default: FMS_RC

FMS_RC:  $(OBJECTS)
	 $(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)

clean:
	rm -f *.o *.mod *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.F.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
