set xlabel 'ω (J)'
set ylabel 'κ (J^2τ)'

plot for [j=30:40:1] 'Cofomega/xsize40ysize40/fieldsize20/cofomegaD'.j.'cJ.dat' u 1:($2/0.3**2) w l t sprintf("D = %d cJ",j)
