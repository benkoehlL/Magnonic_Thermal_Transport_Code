set xlabel 'ω (J)'
set ylabel 'κ (J^2τ)'

plot for [j=60:80:5] 'Cofomega/Tdepscat/T30cJ/xsize50ysize30fieldsize8/cofomegaD'.j.'cJ.dat' u 1:($2*(1-0.5*$1+1/6*$1*$1-$1*$1*$1/24)) w l t sprintf("D = %d cJ",j) 
#:(($2)/($1)/0.3*(1.0-exp(-$1/0.3)))
