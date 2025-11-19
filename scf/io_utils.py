# scf/io_utils.py
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
"""

from pathlib import Path


class IOUtils:
    """
    Clase estática con métodos de entrada/salida de datos.

    Contiene funciones para guardar densidades electrónicas en archivos de
    texto y generar scripts de gnuplot para graficar densidad y
    convergencia SCF.
    """

    @staticmethod
    def save_density(x, rho, filename="rho_vs_x.dat"):
        """
        Guarda la densidad electrónica ρ(x) en un archivo de texto.

        Permite exportar la densidad calculada para graficarla o analizarla posteriormente.

        Parámetros
        ----------
        x : array_like
            Vector de posiciones x de la malla.
        rho : array_like
            Vector de densidad correspondiente a cada posición x.
        filename : str, opcional
            Nombre del archivo de salida. Por defecto 'rho_vs_x.dat'.

        Devuelve
        -------
        None

        Notas
        -----
        Se escribe un archivo de texto con dos columnas: x y ρ(x), separado por espacios.
        Si el archivo ya existía, se sobrescribe.

        Ejemplo
        -------
        >>> IOUtils.save_density([0,1,2], [0.1,0.2,0.3], filename='rho.dat')
        """
        with open(filename, "w", encoding="utf-8") as f:
            for xi, rhoi in zip(x, rho):
                f.write(f"{xi} {rhoi}\n")

    @staticmethod
    def write_gnuplot_density_script(E: float, potential_type: int, output_gp="plot.gp"):
        """
        Genera un script de gnuplot para graficar la densidad ρ(x) según el potencial.

        Parámetros
        ----------
        E : float
            Energía del sistema correspondiente a la densidad.
        potential_type : int
            Tipo de potencial:
            1 → V(x) = α |x|
            2 → V(x) = β (x² - a²)²
            3 → V(x) = 0.5 m ω² x²
            4 → V(x) pozo cuadrado
        output_gp : str, opcional
            Nombre del archivo de salida del script gnuplot. Por defecto 'plot.gp'.

        Devuelve
        -------
        None

        Notas
        -----
        El script generado crea un gráfico EPS con ρ(x) vs x.
        El título del gráfico incluye el tipo de potencial y la energía E.

        Ejemplo
        -------
        >>> IOUtils.write_gnuplot_density_script(E=0.5, potential_type=3)
        """
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
        """
        Genera un script de gnuplot para graficar la convergencia del ciclo SCF.

        Parámetros
        ----------
        output_gp : str, opcional
            Nombre del archivo de salida del script gnuplot. Por defecto 'convergencia.gp'.

        Devuelve
        -------
        None

        Notas
        -----
        El script genera un gráfico EPS de |Δρ| vs iteración en escala logarítmica.
        Permite visualizar el progreso y la estabilidad de la convergencia SCF.

        Ejemplo
        -------
        >>> IOUtils.write_gnuplot_convergence_script()
        """
        with open(output_gp, "w", encoding="utf-8") as g:
            g.write("set term epscairo color\n")
            g.write("set output 'convergencia.eps'\n")
            g.write("set title 'Convergencia SCF'\n")
            g.write("set xlabel 'Iteración'\nset ylabel '|Δρ|'\nset logscale y\n")
            g.write("plot 'convergencia.dat' u 1:2 w lp lw 2 lc rgb '#d62728' t 'error'\n")

