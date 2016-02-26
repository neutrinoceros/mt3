#Plot fitting amplitude between series and Observation

reset
set terminal pdf color
set output "pictures/fit_amplitude_ser_obs.pdf"
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
plot datafile u 5:3 w l , datafile u 5:1 w l

#Plotting dY
set origin 0.0,0.0
set notitle
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set xlabel "t (JC)"
set ylabel "dY (mas)"
plot datafile u 5:4 w l, datafile u 5:2 w l
unset multiplot
