unset key
unset grid
unset xlabel
unset ylabel
unset zlabel
set ytics
set xyplane at 0
set zrange [0 to 1.6]
set yrange [0 to 1]
set xtics (0.5,1,1.5,2)
splot for [i=1:11] "fulldata.out" u 1:2*i:2*i+1 w lines