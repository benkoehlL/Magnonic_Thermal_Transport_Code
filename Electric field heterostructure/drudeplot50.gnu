set xlabel 'E (J)'
set ylabel 'κ (J^2τ)'

T=20

plot for [j=0:50:2] 'Drudeweightoffield/Tdepscat/T'.T.'cJ/drudeweightxsize50ysize30fieldsize'.j.'.dat' u 1:(($2*2)/2/(T*0.01)**2) t sprintf("l  = %d",j)
