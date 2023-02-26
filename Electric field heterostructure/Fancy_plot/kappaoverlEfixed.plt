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
set xrange [0:50] 
set yrange [0:1]
# tics 
set xtics 0,10,50  scale 1.5 ,0.75 format ''
set mxtics 5
set xtics add ('\psb{0}' 0 , '\psb{10}' 10, '\psb{20}' 20, '\psb{30}' 30, '\psb{30}' 30, '\psb{40}' 40, '\psb{50}' 50)
set ytics 0,0.2,1.0 scale 1.5 ,0.75 format '' 
set mytics 2
set ytics add ('\psb{0.2}' 0.2, '\psb{0.4}' 0.4, '\psb{0.6}' 0.6, '\psb{0.8}' 0.8, '\psb{1.0}' 1.0)
# axes-labels 
set xlabel '\psb{$l/a$}'
set ylabel '\psb{$\kappa_{xx}/\kappa_{xx,0}$}' offset 0,0
#labels 
set label '\psb{E = 0.1J}' at graph 0.1,0.91
# key 
set key Left  width -3 at graph 0.9,0.95 samplen 3.5 spacing 1.0
#plot 

E=10

plot for [t=10:30:10] 'Drudeoverl/kappaoverlT'.t.'cJE'.E.'cJ.dat' u 1:2 w l lw 2 t sprintf("T  = %2.1f J",t*0.01)
