#Plot

reset
set terminal pdf color
set output "pictures/fit_amplitude.pdf"
set nokey
set multiplot 
set size 1, 0.5

set origin 0.0,0.5
set title "The Changement of Amplitude X and Y"
set noxlabel 
set ylabel "X"
plot 'amplitude.dat' u 3:1 with linespoints lc 1 #Plot X amplitude

set origin 0.0,0.0
set notitle
set xlabel "Number of Amplitude"
set ylabel "Y"
plot 'amplitude.dat' u 3:2 with linespoints lc 2 #Plot Y amplitude
unset multiplot