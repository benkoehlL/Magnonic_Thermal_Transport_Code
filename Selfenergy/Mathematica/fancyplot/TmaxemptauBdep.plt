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
set yrange [0.2:0.35]
# tics 
set xtics 0,0.1,1  scale 1.5 ,0.75 format ''
set mxtics 5
set xtics add ('\psb{0}' 0 , '\psb{0.2}' 0.2, '\psb{0.4}' 0.4, '\psb{0.6}' 0.6, '\psb{0.8}' 0.8 , '\psb{1.0}' 1.0)
set ytics 0,0.01,0.4 scale 1.5 ,0.75 format '' 
set mytics 2
set ytics add ( '\psb{0.2}' 0.2, '\psb{0.25}' 0.25, '\psb{0.30}' 0.3, '\psb{0.35}' 0.35)
# axes-labels 
set xlabel '\psb{$B/B_S$}'
set ylabel '\psb{$T(\kappa_{max})/J$}' offset 0,0
#labels 
set label '.' at graph -0.19,0.91 textcolor rgb "white"
set label '.' at graph 1.05,1.05 textcolor rgb "white"
# key 
unset key
#plot 

plot 'Bfieldtauempmaximum.dat' u ($2):($1) w l lw 7.5
