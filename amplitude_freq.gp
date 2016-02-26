#Plot the amplitude of dX and dY with frequency

#definition of the environement of the gnuplot
reset
set terminal pdf color
set output "pictures/amplitude_freq.pdf"
set nokey
set multiplot 
set size 1.0, 0.5

#Variable declaration
data = 'amplitude.dat'
#xmin = -10000
#xmax =  10000
#ymin = -0.05
#ymax =  0.04

#Plotting in dX amplitude
set origin 0.0,0.5
#set xrange [xmin:xmax]
#set yrange [ymin:ymax]
set noxlabel
set xlabel "f (rad/JC)"
set ylabel "dX (mas)"
plot data u 5:1 w p pt 7 ps 0.3

#Plotting in dY amplitude
set origin 0.0, 0.0
set notitle
#set xrange [xmin:xmax]
#set yrange [ymin:ymax]
set xlabel "f (rad/JC)"
set ylabel "dY (mas)"
plot data u 5:2 w p pt 7 ps 0.3

unset multiplot
