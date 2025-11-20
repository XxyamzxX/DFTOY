# tests/test_grid.py
# -*- coding: utf-8 -*-
"""
Unit tests para grid.py

POR QUÉ: Grid1D es la base espacial de todo el cálculo. Errores en x o dx
causan integrales incorrectas, normalizaciones erróneas y desalineación
entre potenciales y densidades.
"""

import pytest
import numpy as np
from scf.grid import Grid1D

def test_grid_initialization():
    """
    Verifica que Grid1D inicializa correctamente x y dx.
    
    POR QUÉ: __post_init__ debe crear el array x con np.linspace y calcular dx.
    Si falla, toda función sobre la malla está mal posicionada.
    """
    grid = Grid1D(N=101, x_min=-5.0, x_max=5.0)
    assert grid.x.shape == (101,), "x debe tener N elementos"
    assert grid.x[0] == -5.0, "Primer punto debe ser x_min"
    assert grid.x[-1] == 5.0, "Último punto debe ser x_max"
    assert abs(grid.dx - 0.1) < 1e-9, f"dx incorrecto: {grid.dx}"

def test_grid_uniformity():
    """
    Verifica que los puntos de la malla estén uniformemente espaciados.
    
    POR QUÉ: Si np.linspace falla o se usa mal, los puntos no son equidistantes
    y operadores diferenciales dan resultados incorrectos.
    """
    grid = Grid1D(N=50, x_min=0.0, x_max=1.0)
    diffs = np.diff(grid.x)
    assert np.allclose(diffs, grid.dx, atol=1e-12), "Espaciamiento no uniforme"

def test_grid_small_N():
    """
    Verifica que la malla funcione con N pequeño (caso límite).
    
    POR QUÉ: N=2 es el mínimo para definir un intervalo; verifica robustez
    en casos extremos y fórmula dx = (x_max-x_min)/(N-1).
    """
    grid = Grid1D(N=2, x_min=0.0, x_max=1.0)
    assert len(grid.x) == 2
    assert grid.dx == 1.0  # (1-0)/(2-1) = 1
