set term epscairo color
set output 'rho_vs_x.eps'
set title 'Densidad para V(x)=α|x|, E=1.997401'
set xlabel 'x'
set ylabel 'ρ(x)'
plot 'rho_vs_x.dat' u 1:2 w l lw 2 lc rgb '#1f77b4' t 'ρ(x)'
