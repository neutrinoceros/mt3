#Plot the amplitude of dX and dY with frequency

#definition of the environement of the gnuplot
reset
set terminal pdf color
set output "pictures/fit_amplitude.pdf"
set nokey
set multiplot 
set size 1, 0.5

#Plotting in dX amplitude
set origin 0.0,0.5
set title "The Change of Amplitude dX and dY"
set noxlabel 
set ylabel "dX"
plot 'amplitude.dat' u 3:1 w linespoints lc 1

#Plotting in dY amplitude
set origin 0.0,0.0
set notitle
set xlabel "Number of Amplitude"
set ylabel "dY"
plot 'amplitude.dat' u 3:2 w linespoints lc 2	
unset multiplot
