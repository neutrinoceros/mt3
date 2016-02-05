FC    = gfortran
FLAGS = -g -O0 -fcheck=all -fbackslash
#FLAGS = -g -O2 -Wall -Wextra -fcheck=all -fbackslash

all : lsq clean

lsq : arg_nut.o lsq.f90 
	$(FC) $(FLAGS) $^ -llapack -lblas -o lsq

arg_nut.o : arg_nut.f90
	$(FC) -c $^

run: lsq
	./lsq

###
clean:
	rm *.o *.mod *~ > /dev/null
