OBJ = help.o parse.o extendingBZ.o test.o recip.o \
      plindx.o  vel_exp.o funcs.o ener_read.o \
      kubo_berry_curvature.o kubo_berry_curvature_nonabelian.o \
      get_ext_variable.o get_sorted_xrvari.o \
      write_result.o \
      main.o

# compiler
F90  = mpiifort -fpp -assume byterecl -qmkl -DMPI_USE

 
main : $(OBJ)
	$(F90) $(OBJ) -L/appl/intel2021/oneapi/mkl/2022.0.2/lib/intel64/ -lmpi  -o v2qgt.x

.SUFFIXES: .o .f90

.f90.o :
	$(F90) -c $*.f90

clean :
	rm -f *.o *.mod *~ v2qgt.x

