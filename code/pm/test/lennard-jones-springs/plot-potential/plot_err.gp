
#set terminal jpeg truecolor 
#set output "plot.jpg"

set nokey
set xlabel "Distance (nm)"
set ylabel "Energy (pN.nm)"

set label "x" at 0.3816,-3.59824e-01

#set xr [0.3:0.5]
set yr [0.0:100]

plot "energy.dat" using 1:4 with lines

pause -1


