set grid
#set xlabel 'Temperatur (J)'
set xtics -10,.5,-1
set xtics add ('G' 0 , 'M' 4.44, 'X1' 7.58, 'X2' 12.03, 'G' 15.17)
set ylabel 'tau/J'
#plot for [i=1:7:2] 'interpolateddecayrateplot/interpolateddecayrateplotB'.i.'00mBs.dat' w l t sprintf("B = %d00 mB_S",i),
plot for [j=800:900:250] 'interpolateddecayrateplot/interpolateddecayrateplotB'.j.'mBs.dat' w l lw 7 t sprintf("B = %d mB_S",j)
