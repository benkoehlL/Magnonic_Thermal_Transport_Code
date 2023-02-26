set xlabel 'ω (J)'
set ylabel 'κ (J^2τ)'

plot for [j=35:45:1] 'Cofomega/cofomegaD'.j.'cJxsize50ysize30fieldsize22.dat' u 1:($2/0.3**2) w l t sprintf("D = %d cJ",j)
