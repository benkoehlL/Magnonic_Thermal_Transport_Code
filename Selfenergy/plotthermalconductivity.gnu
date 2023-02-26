set grid
set xlabel 'Temperatur (J)'
set ylabel 'heat conductivity (arb. units)'
set xtics 0,0.1,1.0
#plot for [i=1:7:2] 'Thermalconductivity/thermalconductivityB'.i.'00mBs.dat' w l t sprintf("B = %d00 mB_S",i),
#plot for [j=800:975:100] 'Thermalconductivity/thermalconductivityB'.j.'mBs.dat' w p t sprintf("B = %d mB_S",j)
plot 'SelfenergyB90cBsvariousTemperatures/thermalconductivityB90cBs.dat' w l lw 5 
