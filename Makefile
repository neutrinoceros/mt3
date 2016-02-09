FC    = gfortran
FLAGS = -g -O0 -fcheck=all -fbackslash
LFLAG = -llapack -lblas 
#FLAGS = -g -O2 -Wall -Wextra -fcheck=all -fbackslash

all : lsq 

lsq : mod_matrix.o arg_nut.o lsq.f90 
	$(FC) $(FLAGS) $^ -o lsq $(LFLAG)

arg_nut.o : arg_nut.f90
	$(FC) -c $^

mod_matrix.o : mod_matrix.f90
	$(FC) -c $^ 

run: lsq
	./lsq

###
clean:
	rm *.o *.mod *~ > /dev/null
