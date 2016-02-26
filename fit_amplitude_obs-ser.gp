#Plot the different between fitting amplitude in observation and series

reset
set terminal pdf color
set output "pictures/fit_amplitude_obs-ser.pdf"
set nokey
set multiplot
set size 1, 0.5

#Variable declaration
datafile='series.dat'
xmin = -0.21
xmax =  0.16
ymin = -140
ymax =  60

#Plotting dX
set origin 0.0,0.5
set notitle
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set noxlabel
set ylabel "dX (mas)"
plot datafile u 5:($3-$1) w l

#Plotting dY
set origin 0.0,0.0
set notitle
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set xlabel "time (JC)"
set ylabel "dY (mas)"
plot datafile u 5:($4-$2) w l
unset multiplot
