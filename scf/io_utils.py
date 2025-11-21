# -*- coding: utf-8 -*-
"""
io_utils
========

Proporciona utilidades de entrada/salida para DFToy, incluyendo funciones
para guardar densidades electrónicas y generar scripts de gnuplot para
visualizar densidades y convergencia SCF.

Ejemplo
-------
>>> from scf.io_utils import IOUtils
>>> IOUtils.save_density([0,1,2], [0.1,0.2,0.3], filename='rho.dat')
>>> IOUtils.write_gnuplot_density_script(E=0.5, potential_type=3)
>>> IOUtils.write_gnuplot_convergence_script()
>>> IOUtils.run_gnuplot('plot.gp', 'convergencia.gp')
"""

from pathlib import Path
import subprocess

class IOUtils:
    """
    Clase estática con métodos de entrada/salida de datos.

    Contiene funciones para guardar densidades electrónicas en archivos de
    texto y generar scripts de gnuplot para graficar densidad y
    convergencia SCF.
    """

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

    @staticmethod
    def open_eps(filename: str):
        """
        Abre un archivo EPS usando el visualizador predeterminado del sistema.
        """
        import subprocess, platform
        if platform.system() == 'Darwin':       # macOS
            subprocess.run(['open', filename])
        elif platform.system() == 'Windows':    # Windows
            subprocess.run(['start', filename], shell=True)
        else:                                   # Linux
            subprocess.run(['xdg-open', filename])

    @staticmethod
    def run_gnuplot(*scripts):
        """
        Ejecuta uno o varios scripts de gnuplot.

        Parámetros
        ----------
        scripts : str
            Nombres de archivos de script de gnuplot a ejecutar.

        Ejemplo
        -------
        >>> IOUtils.run_gnuplot('plot.gp', 'convergencia.gp')
        """
        for script in scripts:
            subprocess.run(["gnuplot", script])
