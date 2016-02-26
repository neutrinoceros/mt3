#Plot the amplitude of dX and dY with frequency

#definition of the environement of the gnuplot
reset
set terminal pdf color
set output "pictures/amplitude_freq.pdf"
set nokey
set multiplot title "The Change of Amplitude dX and dY"
set size 0.5, 0.9

#Variable declaration
data = 'amplitude.dat'
xmin = -10000
xmax =  10000
ymin = -0.05
ymax =  0.04

#Plotting in dX amplitude
set origin 0.5, 0.02
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set noxlabel
set xlabel "Frequency (rad/JC)"
set ylabel "dX (mas)"
plot data u 5:1:3 w errorb pt 7 ps 0.1

#Plotting in dY amplitude
set origin 0.0, 0.02
set notitle
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set xlabel "Frequency (rad/JC)"
set ylabel "dY (mas)"
plot data u 5:2:4 w errorb pt 7 ps 0.1

unset multiplot
