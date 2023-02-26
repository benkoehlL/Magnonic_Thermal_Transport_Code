set terminal epslatex size 4.3,3.0 color colortext  
set out 'plot.tex'
#
# plotting style
#
set multiplot
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
set yrange [0:1.35]
# tics 
set xtics 0,0.25,3.1  scale 1.5 ,0.75 format ''
set mxtics 5
set xtics add ('\psb{0}' 0 , '\psb{0.5}' 0.5, '\psb{1}' 1, '\psb{1.5}' 1.5, '\psb{2}' 2, '\psb{2.5}' 2.5, '\psb{3}' 3.0)
set ytics 0,0.1,1.5 scale 1.5 ,0.75 format '' 
set mytics 2
set ytics add ('\psb{0.2}' 0.2, '\psb{0.4}' 0.4, '\psb{0.6}' 0.6, '\psb{0.8}' 0.8, '\psb{1.0}' 1.0, '\psb{1.2}' 1.2)
# axes-labels 
set xlabel '\psb{$T/J$}'
set ylabel '\psb{$\kappa/J$}' offset 0,0
#labels 
set label '.' at graph -0.19,0.91 textcolor rgb "white"
set label '.' at graph 1.05,1.05 textcolor rgb "white"
# key 
set key Left  width -3 at graph 0.975,0.375 samplen 3.5 spacing 1.0
#plot 

plot 'BfieldtauempTdep.dat' u 1:($2) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.0), 'BfieldtauempTdep.dat' u 1:($5) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.3), 'BfieldtauempTdep.dat' u 1:($7) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.5), 'BfieldtauempTdep.dat' u 1:($10) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.8)

#inset
unset label

unset tics

set size 0.45,0.45
set origin 0.5,0.5
set xrange [0:0.10]
set yrange [0:0.35]
set xtics 0,0.025,0.3  scale 1.5 ,0.75 format ''
set mxtics 2
set xtics add ('\psa{0}' 0 , '\psa{0.05}' 0.05, '\psa{0.1}' 0.1)
set ytics 0,0.05,0.5 scale 1.5 ,0.75 format '' 
set mytics 2
set ytics add ('\psa{0.1}' 0.1, '\psa{0.2}' 0.2, '\psa{0.3}' 0.3)
set xlabel ''
set ylabel '' offset 0,0
unset key

plot 'BfieldtauempTdep.dat' u 1:(0.5*$2) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.0), 'BfieldtauempTdep.dat' u 1:($5) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.3), 'BfieldtauempTdep.dat' u 1:($7) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.5), 'BfieldtauempTdep.dat' u 1:($10) w l lw 7.5 t sprintf("$B$ = %0.1f $B_S$",0.8)

unset multiplot
