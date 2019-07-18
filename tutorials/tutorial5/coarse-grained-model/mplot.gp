#set terminal jpeg truecolor
#set output "plot.jpg"

set key outside
set xlabel "Time"
set ylabel "Energy (pN.nm)"

set style line 1 lt 1 lc 1 lw 2
set style line 2 lt 1 lc 2 lw 2
set style line 3 lt 1 lc 3 lw 2

plot "m1rsim_ljspring_energy.dat" using 1:2 with lines ls 1 title "c-alpha", \
     "m2rsim_ljspring_energy.dat" using 1:2 with lines ls 2 title "LJ-sc", \
     "m3rsim_ljspring_energy.dat" using 1:2 with lines ls 3 title "LJ-all" \

pause -1


