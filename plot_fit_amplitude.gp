#Plot

reset
set multiplot 
set size 1, 0.5

#set origin 0.0,0.5
#set title "The Changement of Amplitude X and Y"
#set noxlabel 
#set ylabel "X"
#plot 'amplitude.txt' u 3:1 with linespoints lc 1 #Plot X amplitude

#set origin 0.0,0.0
#set notitle
#set xlabel "Number of Amplitude"
#set ylabel "Y"
#plot 'amplitude.txt' u 3:2 with linespoints lc 2 #Plot Y amplitude

file="series.txt"
set origin 0.0,0.5
set title "The comparison of dX and dY between model and observation"
set noxlabel 
set ylabel "dX"
plot file u 5:1 w l, file u 5:3 w l 

set origin 0.0,0.0
set notitle
set xlabel "Time"
set ylabel "dY"
plot file u 5:2  w l, file u 5:4 w l  

unset multiplot
