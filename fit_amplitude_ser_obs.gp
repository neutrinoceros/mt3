#Plot of series versus observation

set terminal pdf color
set output "pictures/fit_amplitude_ser_obs.pdf"
set multiplot
set size 1, 0.4

#Variable declaration
datafile='series.dat'
xmin = -0.05
xmax =  0.16
ymin = -1
ymax =  1

#Plotting dX
set origin 0.0,0.55
set notitle
set key at 0.163,-4.8
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set noxlabel
set ylabel "dX (mas)"
plot datafile u 5:3 w l title "série ajustée", datafile u 5:1 w l title "données observationnelle"

#Plotting dY
set origin 0.0,0.15
set notitle
set nokey
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set xlabel "t (JC)"
set ylabel "dY (mas)"
plot datafile u 5:4 w l title "série ajustée", datafile u 5:2 w l title "données observationnelle"

unset multiplot
