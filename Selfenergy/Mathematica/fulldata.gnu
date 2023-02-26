unset key
unset grid
unset xlabel
unset ylabel
unset zlabel
set xyplane at 0
set zrange [0 to 1.5]
set ytics (0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
set xtics (0.25,0.5,0.75,1)
splot for [i=1:110:10] "fulldata.out" u 1:2*i:2*i+1 w lines