FC    = gfortran
FLAGS = -g -O0 -fcheck=all -fbackslash
LFLAG = -llapack -lblas 
#FLAGS = -g -O2 -Wall -Wextra -fcheck=all -fbackslash

all : lsq lsq_450

lsq : mod_matrix.o arg_nut.o lsq.f90 
	$(FC) $(FLAGS) $^ -o lsq $(LFLAG)

lsq_450 : mod_matrix.o arg_nut.o lsq_450_day_term.f90
	$(FC) $(FLAGS) $^ -o lsq_450 $(LFLAG)

arg_nut.o : arg_nut.f90
	$(FC) -c $^

mod_matrix.o : mod_matrix.f90
	$(FC) -c $^ 

run: lsq
	./lsq

###
clean:
	rm *.o *.mod *~ > /dev/null
