unset key
unset grid
unset xlabel
unset ylabel
unset zlabel
set xyplane at 0
set zrange [0 to 1.5]
set ytics (0.,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1)
set xtics (1,2,3)
splot for [i=1:20:1] "smallbdata.out" every 5 u 1:2*i:2*i+1