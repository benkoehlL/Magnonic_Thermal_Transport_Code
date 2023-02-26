set terminal epslatex size 4.3,3.0 color colortext  
set out 'plot.tex'
#
# plotting style
#
set style line 1 lc rgb "dark-blue"  dt 1 lw 1.5 pt 7
set style line 2 lc rgb "dark-violet" dt 2 lw 3.0 pt 5
#
# configure plot 
#
unset title  
set border lw 2.0
# margins 
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

# Plotting
# 
#ranges 
set xrange [0:1] 
set yrange [0.35:1.35]
# tics 
set xtics 0,0.1,1  scale 1.5 ,0.75 format ''
set mxtics 5
set xtics add ('\psb{0}' 0 , '\psb{0.2}' 0.2, '\psb{0.4}' 0.4, '\psb{0.6}' 0.6, '\psb{0.8}' 0.8 , '\psb{1.0}' 1.0)
set ytics 0,0.1,1.4 scale 1.5 ,0.75 format '' 
set mytics 2
set ytics add ( '\psb{0.4}' .4, '\psb{0.6}' .6, '\psb{0.8}' 0.8, '\psb{1.0}' 1.0, '\psb{1.2}' 1.2 )
# axes-labels 
set xlabel '\psb{$B/B_S$}'
set ylabel '\psb{$\kappa(B)/\kappa(0)$}' offset 0,0
#labels 
set label '.' at graph -0.19,0.91 textcolor rgb "white"
set label '.' at graph 1.05,1.05 textcolor rgb "white"
# key 
set key Left  width -3 at graph 0.375,0.95 samplen 3.5 spacing 1.0
#plot 

plot 'kappafulltauconstBdepnormed.dat' u 1:($3) w l lw 7.5 t sprintf("$T$ = %0.2f $J$",0.1), 'kappafulltauconstBdepnormed.dat' u 1:($5) w l lw 7.5 t sprintf("$T$ = %0.2f $J$",0.2), 'kappafulltauconstBdepnormed.dat' u 1:($7) w l lw 7.5 t sprintf("$T$ = %0.2f $J$",0.3), 'kappafulltauconstBdepnormed.dat' u 1:($9) w l lw 7.5 t sprintf("$T$ = %0.2f $J$",0.4), 'kappafulltauconstBdepnormed.dat' u 1:($11) w l lw 7.5 t sprintf("$T$ = %0.2f $J$",0.5)
