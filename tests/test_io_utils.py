# tests/test_io_utils.py
# -*- coding: utf-8 -*-
"""
Unit tests para io_utils.py

POR QUÉ: IOUtils guarda resultados a disco; errores causan pérdida de datos,
archivos corruptos o gnuplot scripts inválidos.
"""

import pytest
import os
from pathlib import Path
from scf.io_utils import IOUtils

def test_save_density_creates_file(tmp_path):
    """
    Verifica que save_density cree archivo correctamente.
    
    POR QUÉ: Si el archivo no se crea o tiene formato incorrecto,
    la GUI no puede leerlo y gnuplot falla.
    """
    filename = tmp_path / "test_rho.dat"
    x = [0.0, 1.0, 2.0]
    rho = [0.5, 1.0, 0.5]
    IOUtils.save_density(x, rho, filename=str(filename))
    
    assert filename.exists(), "Archivo no fue creado"
    lines = filename.read_text().strip().split('\n')
    assert len(lines) == 3, f"Debe tener 3 líneas, tiene {len(lines)}"

def test_save_density_format(tmp_path):
    """
    Verifica formato correcto: "x rho" en cada línea.
    
    POR QUÉ: Gnuplot y GUI esperan dos columnas separadas por espacio;
    formato incorrecto causa errores al parsear.
    """
    filename = tmp_path / "test_rho.dat"
    x = [1.0, 2.0]
    rho = [0.3, 0.7]
    IOUtils.save_density(x, rho, filename=str(filename))
    
    with open(filename, 'r') as f:
        line = f.readline()
        parts = line.split()
        assert len(parts) == 2, "Cada línea debe tener 2 columnas"
        assert float(parts[0]) == 1.0, "Primera columna debe ser x"
        assert float(parts[1]) == 0.3, "Segunda columna debe ser rho"

def test_write_gnuplot_density_script(tmp_path):
    """
    Verifica que se cree script de gnuplot para densidad.
    
    POR QUÉ: El script debe contener comandos válidos de gnuplot;
    errores de sintaxis causan que gnuplot no grafique.
    """
    filename = tmp_path / "plot.gp"
    IOUtils.write_gnuplot_density_script(E=1.5, potential_type=3, output_gp=str(filename))
    
    assert filename.exists(), "Script gnuplot no fue creado"
    content = filename.read_text()
    assert "set term" in content, "Falta comando 'set term'"
    assert "set output" in content, "Falta comando 'set output'"
    assert "plot" in content, "Falta comando 'plot'"
    assert "E=1.5" in content or "E=1.500000" in content, "Falta energía en título"

def test_write_gnuplot_convergence_script(tmp_path):
    """
    Verifica script de gnuplot para convergencia.
    
    POR QUÉ: El script debe incluir 'set logscale y' para graficar
    error en escala logarítmica; sin esto, convergencia no es visible.
    """
    filename = tmp_path / "conv.gp"
    IOUtils.write_gnuplot_convergence_script(output_gp=str(filename))
    
    assert filename.exists(), "Script convergencia no fue creado"
    content = filename.read_text()
    assert "logscale" in content, "Falta escala logarítmica"
    assert "convergencia.dat" in content, "Falta referencia a archivo de datos"
