set grid
set xlabel 'B (B_S)'
set ylabel 'thermal conductivity (J)'
plot for [i=10:30:10] 'Bfieldthermalconductivity/bfieldthermalconductivityT'.i.'cJ.dat' u 1:2 w l t sprintf("T = %d cJ",i), for [i=10:30:10] 'Bfieldthermalconductivity/bfieldthermalconductivityT'.i.'cJ.dat' u 1:3 w l t sprintf("T = %d cJ",i)
