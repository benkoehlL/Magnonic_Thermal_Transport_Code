set xlabel 'E (J)'
set ylabel 'κ (J^2τ)'

T=10

plot for [j=2:40:2] 'Drudeweightoffield/T'.T.'cJ/drudeweightxsize40ysize40fieldsize'.j.'.dat' u 1:(($2*2)/2/(T*0.01)**2) t sprintf("l  = %d",j), 'Drudeweightoffield/T'.T.'cJ/drudeanalyticxsize40ysize40.dat' u 1:(($2*2)/2/(T*0.01)**2) t sprintf("l  = 40")

