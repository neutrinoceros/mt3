#Plot fitting amplitude between series - MHB and Observation

reset
set terminal pdf color
set output "pictures/fit_amplitude_ser_obs_mhb.pdf"
set nokey
set multiplot
set size 1, 0.5

#Plotting dX
set origin 0.0,0.5
set title "Comparision between the observation and the series"
set noxlabel
set ylabel "dX"
plot 'series.txt' u 5:1 w l , 'series.txt' u 5:3 w l

#Plotting dY
set origin 0.0,0.0
set notitle
set xlabel "time"
set ylabel "dY"
plot 'series.txt' u 5:2 w l, 'series.txt' u 5:4 w l
unset multiplot
