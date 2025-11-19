set term epscairo color
set output 'rho_vs_x.eps'
set title 'Densidad para V(x)=0.5 m ω^2 x^2, E=3.289070'
set xlabel 'x'
set ylabel 'ρ(x)'
plot 'rho_vs_x.dat' u 1:2 w l lw 2 lc rgb '#1f77b4' t 'ρ(x)'
