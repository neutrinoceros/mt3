FC    = gfortran
FLAGS = -g -O0 -fcheck=all -fbackslash
LFLAG = -llapack -lblas 
#FLAGS = -g -O2 -Wall -Wextra -fcheck=all -fbackslash

all : lsq_457 lsq 

lsq : mod_matrix.o arg_nut.o mod_series.o lsq.f90 
	$(FC) $(FLAGS) $^ -o lsq.exe $(LFLAG)

lsq_457 : arg_nut.o mod_matrix.o mod_lsq.o mod_series.o lsq_457_day_term.f90
	$(FC) $(FLAGS) $^ -o lsq_457.exe $(LFLAG)

arg_nut.o : arg_nut.f90
	$(FC) -c $^

mod_matrix.o : mod_matrix.f90
	$(FC) -c $^ 

mod_series.o : mod_series.f90
	$(FC) -c $^

mod_lsq.o : mod_lsq.f90
	$(FC) -c $^

run: lsq_457 lsq
	./lsq_457.exe
	./lsq.exe

###
clean :
	rm *.dat *.o *.mod *.exe 2> /dev/null

rebuild : clean all


rerun : clean run

