FC    = gfortran
FLAGS = -g -O0
#FLAGS = -g -O2 -Wall -Wextra -fcheck=all -fbackslash

all: lapack

lapack: dgetri.o dtrti2.o dtrtri.o ieeeck.o ilaenv.o iparmq.o lsame.o xerbla.o

dgetri.o : lapack/double/* lapack/lapack_routine/*
	$(FC) $(FLAGS) -c lapack/double/dgetri.f

dtrti2.o : lapack/double/dtrti2.f
	$(FC) $(FLAGS) -c $^

dtrtri.o : lapack/double/dtrtri.f
	$(FC) $(FLAGS) -c $^

ieeeck.o : lapack/lapack_routine/ieeeck.f
	$(FC) $(FLAGS) -c $^

ilaenv.o : lapack/lapack_routine/ilaenv.f
	$(FC) $(FLAGS) -c $^

iparmq.o : lapack/lapack_routine/iparmq.f
	$(FC) $(FLAGS) -c $^

lsame.o : lapack/lapack_routine/lsame.f
	$(FC) $(FLAGS) -c $^

xerbla.o : lapack/lapack_routine/xerbla.f
	$(FC) $(FLAGS) -c $^

lsq : arg_nut.o lsq.f90 
	$(FC) $(FLAGS) $^ -llapack -lblas -o lsq

arg_nut.o : arg_nut.f90
	$(FC) -c $^
run: lsq
	./lsq

###
clean:
	rm *.o *.mod *~ > /dev/null
