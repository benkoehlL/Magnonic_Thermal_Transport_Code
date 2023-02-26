set xlabel 'E (J)'
set ylabel 'κ (J^2τ)'

T=30

plot for [j=2:20:2] 'Drudeweightoffield/T'.T.'cJ/drudeweightxsize20ysize20fieldsize'.j.'.dat' u 1:(($2*2)/2/(T*0.01)**2) t sprintf("l  = %d",j)
