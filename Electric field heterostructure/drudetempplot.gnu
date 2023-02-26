set xlabel 'T (J)'
set ylabel 'κ (1/J)'

f=26

plot 'DrudeweightofT/Tdepscat/drudeweightnumD1cJxsize50ysize30fieldsize0.dat' u 1:($2/$1**2) w l t sprintf("D = 0 cJ"), for [j=10:20:10] 'DrudeweightofT/Tdepscat/drudeweightnumD'.j.'cJxsize50ysize30fieldsize'.f.'.dat' u 1:($2/$1**2) w l t sprintf("D = %d cJ, l=26",j), for [j=10:20:10] 'DrudeweightofT/Tdepscat/drudeweightnumD'.j.'cJxsize50ysize30fieldsize24.dat' u 1:($2/$1**2) w l t sprintf("D = %d cJ, l=24",j), for [j=10:20:10] 'DrudeweightofT/Tdepscat/drudeweightnumD'.j.'cJxsize50ysize30fieldsize50.dat' u 1:($2/$1**2) w l t sprintf("D = %d cJ, l=50",j), for [j=10:20:10] 'DrudeweightofT/Tdepscat/drudeweightnumD'.j.'cJxsize50ysize30fieldsize2.dat' u 1:($2/$1**2) w l t sprintf("D = %d cJ, l=2",j)

