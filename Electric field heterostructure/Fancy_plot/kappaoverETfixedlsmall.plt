T=30
N=(4.04555*2)/2/(T*0.01)**2


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
set xrange [0:1] 
set yrange [0:1]
# tics 
set xtics 0,0.2,1.0  scale 1.5 ,0.75 format ''
set mxtics 2
set xtics add ('\psb{0}' 0 , '\psb{0.2}' 0.2, '\psb{0.4}' 0.4, '\psb{0.6}' 0.6, '\psb{0.8}' 0.8, '\psb{1.0}' 1.0)
set ytics 0,0.2,1.0 scale 1.5 ,0.75 format '' 
set mytics 2
set ytics add ('\psb{0.2}' 0.2, '\psb{0.4}' 0.4, '\psb{0.6}' 0.6, '\psb{0.8}' 0.8, '\psb{1.0}' 1.0)
# axes-labels 
set xlabel '\psb{$E/J$}'
set ylabel '\psb{$\kappa_{xx}/\kappa_{xx,0}$}' offset 0,0
#labels 
set label '\psb{T = 0.3J}' at graph 0.35,0.9
# key 
set key Left  width -3 at graph 0.9,0.92 samplen 3.5 spacing 1.0
#plot 

plot 'Drudeweightoffield/Tdepscat/T'.T.'cJ/drudeweightxsize50ysize30fieldsize2.dat' u 1:(($2*2)/2/(T*0.01)**2/N) w l lw 2 t sprintf("l  = 2a"), 'Drudeweightoffield/Tdepscat/T'.T.'cJ/drudeweightxsize50ysize30fieldsize4.dat' u 1:(($2*2)/2/(T*0.01)**2/N) w l lw 2 t sprintf("l  = 4a"), 'Drudeweightoffield/Tdepscat/T'.T.'cJ/drudeweightxsize50ysize30fieldsize8.dat' u 1:(($2*2)/2/(T*0.01)**2/N) w l lw 2 t sprintf("l  = 8a"), 'Drudeweightoffield/Tdepscat/T'.T.'cJ/drudeweightxsize50ysize30fieldsize50.dat' u 1:(($2*2)/2/(T*0.01)**2/N) w l lw 2 t sprintf("l  = 50a")
