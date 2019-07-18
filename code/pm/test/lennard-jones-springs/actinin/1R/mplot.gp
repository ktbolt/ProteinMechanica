#set terminal jpeg truecolor
#set output "plot.jpg"

set key outside
set xlabel "Time"
set ylabel "Energy (pN.nm)"

set style line 1 lt 1 lc 1 lw 2
set style line 2 lt 1 lc 2 lw 2
set style line 3 lt 1 lc 3 lw 2
set style line 4 lt 1 lc 4 lw 2

plot "rsim1_energy.dat" using 1:2 with lines ls 1 title "model1", \
     "rsim2_energy.dat" using 1:2 with lines ls 2 title "model2", \
     "rsim3_energy.dat" using 1:2 with lines ls 3 title "model3", \
     "rsim4_energy.dat" using 1:2 with lines ls 4 title "model4"

pause -1


