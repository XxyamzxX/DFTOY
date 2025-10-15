# -*- coding: utf-8 -*-
from pathlib import Path
import subprocess

class IOUtils:
    @staticmethod
    def save_density(x, rho, filename="rho_vs_x.dat"):
        with open(filename, "w") as f:
            for xi, rhoi in zip(x, rho):
                f.write(f"{xi} {rhoi}\n")

    @staticmethod
    def open_eps(eps_filename: str):
        p = Path(eps_filename)
        if not p.exists():
            return
        # Preferir xdg-open; fallback a evince
        try:
            subprocess.Popen(["xdg-open", eps_filename])
        except Exception:
            try:
                subprocess.Popen(["evince", eps_filename])
            except Exception:
                pass

    @staticmethod
    def write_gnuplot_density_script(E: float, potential_type: int, output_gp="plot.gp"):
        title_map = {
            1: "Densidad electronica para V(x) = {/Symbol a}|x|",
            2: "Densidad electronica para V(x) = {/Symbol b}(x^2 - a^2)^2",
            3: "Densidad electronica para V(x) = 0.5 m {/Symbol w}^2 x^2",
            4: "Densidad electronica para V(x) = -V0 (|x|<L/2)",
        }
        title = title_map.get(potential_type, "Densidad electronica")
        with open(output_gp, "w") as gp:
            gp.write("set terminal postscript eps enhanced color font 'Times,18'\n")
            gp.write("set output 'rho_vs_x.eps'\n")
            gp.write("set xlabel 'x'\n")
            gp.write("set ylabel '{/Symbol r}(x)'\n")
            gp.write("set grid\n")
            gp.write(f"set title '{title}'\n")
            gp.write("set style textbox opaque border\n")
            gp.write(f"set label 1 sprintf('E_0 = %.6f J', {E}) at graph 0.05,0.95 left front boxed\n")
            gp.write("plot 'rho_vs_x.dat' with lines lw 2 title '{/Symbol r}(x)'\n")

    @staticmethod
    def write_gnuplot_convergence_script(output_gp="convergencia.gp"):
        with open(output_gp, "w") as gp2:
            gp2.write("set terminal postscript eps enhanced color font 'Times,18'\n")
            gp2.write("set output 'convergencia.eps'\n")
            gp2.write("set title 'Convergencia del ciclo SCF'\n")
            gp2.write("set xlabel 'Iteracion'\n")
            gp2.write("set ylabel '||{/Symbol r} - {/Symbol r}_{old}||'\n")
            gp2.write("set grid\n")
            gp2.write("set logscale y\n")
            gp2.write("set format y '10^{%T}'\n")
            gp2.write("plot 'convergencia.dat' using 1:2 with lines lw 2 lc rgb 'blue' title 'Error de la densidad'\n")

    @staticmethod
    def run_gnuplot(plot_gp="plot.gp", conv_gp="convergencia.gp"):
        try:
            subprocess.run(["gnuplot", plot_gp], check=True)
            subprocess.run(["gnuplot", conv_gp], check=True)
        except Exception:
            print("No se pudo ejecutar gnuplot automáticamente. Ejecute: gnuplot plot.gp ; gnuplot convergencia.gp")
