# scf/io_utils.py
# -*- coding: utf-8 -*-
from pathlib import Path

class IOUtils:
    @staticmethod
    def save_density(x, rho, filename="rho_vs_x.dat"):
        with open(filename, "w", encoding="utf-8") as f:
            for xi, rhoi in zip(x, rho):
                f.write(f"{xi} {rhoi}\n")

    @staticmethod
    def write_gnuplot_density_script(E: float, potential_type: int, output_gp="plot.gp"):
        title_map = {
            1: "Densidad para V(x)=α|x|",
            2: "Densidad para V(x)=β(x^2-a^2)^2",
            3: "Densidad para V(x)=0.5 m ω^2 x^2",
            4: "Densidad para V(x) pozo cuadrado",
        }
        with open(output_gp, "w", encoding="utf-8") as g:
            g.write("set term epscairo color\n")
            g.write("set output 'rho_vs_x.eps'\n")
            g.write(f"set title '{title_map.get(potential_type, 'Densidad electrónica')}, E={E:.6f}'\n")
            g.write("set xlabel 'x'\nset ylabel 'ρ(x)'\n")
            g.write("plot 'rho_vs_x.dat' u 1:2 w l lw 2 lc rgb '#1f77b4' t 'ρ(x)'\n")

    @staticmethod
    def write_gnuplot_convergence_script(output_gp="convergencia.gp"):
        with open(output_gp, "w", encoding="utf-8") as g:
            g.write("set term epscairo color\n")
            g.write("set output 'convergencia.eps'\n")
            g.write("set title 'Convergencia SCF'\n")
            g.write("set xlabel 'Iteración'\nset ylabel '|Δρ|'\nset logscale y\n")
            g.write("plot 'convergencia.dat' u 1:2 w lp lw 2 lc rgb '#d62728' t 'error'\n")
