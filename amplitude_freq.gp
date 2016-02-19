#Plot the amplitude of dX and dY with frequency

#definition of the environement of the gnuplot
reset
set terminal pdf color
set output "pictures/amplitude_freq.pdf"
set nokey
set multiplot 
set size 1, 0.5

#Plotting in dX amplitude
set origin 0.0,0.5
set title "The Change of Amplitude dX and dY"
set noxlabel 
set ylabel "dX"
plot 'amplitude.dat' u 5:1:3 w errorb pt 7 ps 0.7

#Plotting in dY amplitude
set origin 0.0,0.0
set notitle
set xlabel "Frequency"
set ylabel "dY"
plot 'amplitude.dat' u 5:2:4 w errorb pt 7 ps 0.7

unset multiplot
