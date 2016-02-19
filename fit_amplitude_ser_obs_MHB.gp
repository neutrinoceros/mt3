#Plot fitting amplitude between series - MHB and Observation

reset
set terminal pdf color
set output "pictures/fit_amplitude_ser_obs_mhb.pdf"
set nokey
set multiplot
set size 1, 0.5

datafile='series.dat'

#Plotting dX
set origin 0.0,0.5
set title "Comparision between the observation and the series"
set noxlabel
set ylabel "dX (mas)"
plot datafile u 5:1 w l , datafile u 5:3 w l

#Plotting dY
set origin 0.0,0.0
set notitle
set xlabel "time (JC)"
set ylabel "dY (mas)"
plot datafile u 5:2 w l, datafile u 5:4 w l
unset multiplot
