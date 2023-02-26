unset key
unset grid
unset xlabel
unset ylabel
unset zlabel
set xyplane at 0
set zrange [0 to 1.5]
set ytics (0.92,0.94,0.96,0.98,1)
set xtics (1,2,3)
splot for [i=1:20:1] "largebdata.out" every 5 u 1:2*i:2*i+1