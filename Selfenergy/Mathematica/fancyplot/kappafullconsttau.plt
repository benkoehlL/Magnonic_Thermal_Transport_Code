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
set xrange [0:3] 
set yrange [0:0.5]
# tics 
set xtics 0,0.25,3.1  scale 1.5 ,0.75 format ''
set mxtics 5
set xtics add ('\psb{0}' 0 , '\psb{0.5}' 0.5, '\psb{1}' 1, '\psb{1.5}' 1.5, '\psb{2}' 2, '\psb{2.5}' 2.5, '\psb{3}' 3.0)
set ytics 0,0.1,1 scale 1.5 ,0.75 format '' 
set mytics 2
set ytics add ('\psb{0.1}' 0.1, '\psb{0.2}' 0.2, '\psb{0.3}' 0.3, '\psb{0.4}' 0.4, '\psb{0.5}' 0.5, '\psb{0.6}' 0.6, '\psb{0.7}' 0.7, '\psb{0.8}' 0.8, '\psb{0.9}' 0.9)
# axes-labels 
set xlabel '\psb{$T/J$}'
set ylabel '\psb{$\kappa/J^2\tau$}' offset 0,0
#labels 
set label '.' at graph -0.19,0.91 textcolor rgb "white"
set label '.' at graph 1.05,1.05 textcolor rgb "white"
# key 
set key Left  width -3 at graph 0.975,0.325 samplen 3.5 spacing 1.0
#plot 

plot 'BfieldtauconstTdep.dat' u 1:($2) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.0), 'BfieldtauconstTdep.dat' u 1:($5) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.3), 'BfieldtauconstTdep.dat' u 1:($7) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.5), 'BfieldtauconstTdep.dat' u 1:($10) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.8)
