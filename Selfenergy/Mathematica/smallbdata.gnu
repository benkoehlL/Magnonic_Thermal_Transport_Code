unset key
unset grid
unset xlabel
unset ylabel
unset zlabel
set xyplane at 0
set zrange [0 to 20]
set ytics (0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
set xtics (5,10,15)
splot for [i=1:20:1] "smallbdata.out" every 5 u 1:2*i:2*i+1