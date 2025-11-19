set term epscairo color
set output 'convergencia.eps'
set title 'Convergencia SCF'
set xlabel 'Iteración'
set ylabel '|Δρ|'
set logscale y
plot 'convergencia.dat' u 1:2 w lp lw 2 lc rgb '#d62728' t 'error'
