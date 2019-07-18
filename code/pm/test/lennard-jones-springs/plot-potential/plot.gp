
set terminal jpeg truecolor 
set output "plot.jpg"

set nokey
set xlabel "Distance (nm)"
set ylabel "Energy (pN.nm)"

set label "x" at 0.3816,-3.59824e-01

set xzeroaxis
#set xr [0.3:0.5]
#set yr [-0.4:0.2]
set yr [-0.5:0.3]

plot "energy.dat" using 1:2 with lines, "energy.dat" using 1:3 with lines

#pause -1


