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
#
# files to plot 
#
fileA = './DataFiles/data.dat'
#
# functions 
#
fa(x) = (x>0) ? x**2 :1/0  
#
# Plotting
# 
#ranges 
set xrange [0:1.5] 
set yrange [0:58]
# tics 
set xtics 0,0.2,1.6  scale 1.5 ,0.75 format ''
set mxtics 2
set xtics add ('\psb{0}' 0 , '\psb{0.2}' 0.2, '\psb{0.4}' 0.4, '\psb{0.6}' 0.6, '\psb{0.8}' 0.8, '\psb{1}' 1,  '\psb{1.2}' 1.2, '\psb{1.4}' 1.4)
set ytics 0,10,60 scale 1.5 ,0.75 format '' 
set mytics 2
set ytics add ('\psb{10}' 10, '\psb{20}' 20,'\psb{30}' 30,'\psb{40}' 40, '\psb{50}' 50)
# axes-labels 
set xlabel '\psb{$T/J$}'
set ylabel '\psb{$\kappa_{xx}/J$}' offset 0,0
#labels 
set label '\psb{l=50}' at graph 0.3,0.875
# key 
set key at graph 0.9,0.92 samplen 3.5 spacing 1.5
#plot 
f=50

plot 'DrudeweightofT/Tdepscat/drudeweightnumD1cJxsize50ysize30fieldsize0.dat' u 1:($2/$1**2) w l lw 2 t sprintf("D = 0.0 J",0.0), for [j=0:40:10] 'DrudeweightofT/Tdepscat/drudeweightnumD'.j.'cJxsize50ysize30fieldsize'.f.'.dat' u 1:($2/$1**2) w l lw 2 t sprintf("D = %2.1f J",j*0.01)
