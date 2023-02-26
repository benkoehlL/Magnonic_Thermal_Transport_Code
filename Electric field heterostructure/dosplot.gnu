set xlabel 'ω (J)'
set ylabel 'dos (1/J)'

T=20
f=6

plot for [j=10:20:1] 'DOS/xsize50ysize30/T'.T.'cJ/fieldsize'.f.'/D'.j.'cJmod.dat' w l t sprintf("D = %d cJ",j),
#plot for [j=14:50:18] 'DOS/xsize50ysize30/fieldsize'.j.'/D75cJ.dat' w l t sprintf("D = 50, l = %d sites",j)
