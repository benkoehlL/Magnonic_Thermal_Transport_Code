set xlabel 'ω (J)'
set ylabel 'C (1/J)'

T=30
f=4

plot for [j=80:90:1] 'Cofomega/Tdepscat/T'.T.'cJ/xsize50ysize30fieldsize'.f.'/cofomegaD'.j.'cJ.dat' u 1:2 w l t sprintf("D = %d cJ",j)
