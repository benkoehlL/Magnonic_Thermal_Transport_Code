set terminal postscript enhanced
set encoding utf8
set xlabel 'k_x' font ",34"
set ylabel 'k_y' font ",34"
set cblabel '𝞃_{smd} / SJ' font "sfti3400,34"
unset terminal
set tics font ",30"
set xtics offset 0,-1
set ytics offset -1,0
set xlabel offset 0,-4
set ylabel offset -8,0
set cblabel offset 8,0
set tmargin at screen 0.95
set bmargin at screen 0.25
set rmargin at screen 0.77
set lmargin at screen 0.2

plot 'Densityplot/densityplotB850mBs.dat' u 1:2:3 w image
