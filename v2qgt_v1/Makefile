OBJ = help.o parse.o extendingBZ.o test.o recip.o \
      plindx.o  vel_exp.o funcs.o ener_read.o \
      kubo_berry_curvature.o get_ext_variable.o get_sorted_xrvari.o \
      write_result.o \
      main.o

# compiler
F90  = ifort -fpp -assume byterecl -mkl 

 
main : $(OBJ)
	$(F90) $(OBJ)  -o v2geometry.x 
	## cp -f wt.x ../bin

.SUFFIXES: .o .f90

.f90.o :
	$(F90) -c $*.f90

clean :
	rm -f *.o *.mod *~ v2geometry.x

