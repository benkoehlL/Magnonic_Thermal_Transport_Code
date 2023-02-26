set xlabel 'l (sites)'
set ylabel 'κ (J^2τ)'

E=10

plot for [t=10:30:10] 'Drudeoverl/kappaoverlT'.t.'cJE'.E.'cJ.dat' u 1:2 w l t sprintf("T  = %2.2f J",t*0.01)
