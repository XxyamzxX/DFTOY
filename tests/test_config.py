# tests/test_config.py
# -*- coding: utf-8 -*-
"""
Unit tests para config.py

POR QUÉ: Config centraliza constantes físicas y numéricas; errores aquí se propagan
a todo el programa. Estos tests verifican:
1. Valores por defecto razonables (evita inicializaciones accidentales con cero o NaN)
2. Cálculo correcto de dx (espaciamiento de malla)
3. Coeficiente cinético con signo y magnitud correctos (el laplaciano debe ser negativo)
"""

import pytest
from scf.config import Config

def test_config_default_values():
    """
    Verifica que los valores por defecto de Config sean físicamente razonables.
    
    POR QUÉ: Si N, hbar o tolerances son cero, el programa falla silenciosamente
    o produce resultados sin sentido. Este test atrapa errores de inicialización.
    """
    assert Config.N > 0, "N debe ser positivo para definir una malla válida"
    assert Config.x_min < Config.x_max, "x_min debe ser menor que x_max"
    assert Config.tolerance > 0, "tolerance debe ser positiva para convergencia"
    assert Config.max_iter > 0, "max_iter debe ser positivo"
    assert Config.hbar > 0, "hbar debe ser positivo (constante física)"
    assert 0 < Config.mix_alpha < 1, "mix_alpha debe estar en (0,1) para mezcla estable"

def test_dx_calculation():
    """
    Verifica el cálculo de espaciamiento dx = (x_max - x_min)/(N-1).
    
    POR QUÉ: dx es crucial para integrales numéricas y derivadas discretas.
    Un error aquí afecta energías, normalizaciones y convergencia.
    """
    dx = Config.dx(N=201, x_min=-5.0, x_max=5.0)
    expected = 10.0 / 200.0  # (5 - (-5)) / (201 - 1)
    assert abs(dx - expected) < 1e-10, f"dx incorrecto: {dx} vs esperado {expected}"

def test_kinetic_coeff_sign():
    """
    Verifica que el coeficiente cinético sea negativo.
    
    POR QUÉ: El término cinético T = -ℏ²/(2m·dx²) debe ser negativo por definición
    del laplaciano discreto. Un signo erróneo invierte la física del problema.
    """
    coeff = Config.kinetic_coeff(hbar=1.0, m=1.0, dx=0.1)
    assert coeff < 0, "Coeficiente cinético debe ser negativo (laplaciano)"

def test_kinetic_coeff_magnitude():
    """
    Verifica la magnitud del coeficiente cinético.
    
    POR QUÉ: Con dx pequeño, el coeficiente crece mucho; con dx grande, disminuye.
    Esto asegura la fórmula -ℏ²/(2m·dx²) esté implementada correctamente.
    """
    coeff = Config.kinetic_coeff(hbar=1.0, m=2.0, dx=0.5)
    expected = -1.0 / (2.0 * 2.0 * 0.5**2)  # -2.0
    assert abs(coeff - expected) < 1e-10, f"Magnitud incorrecta: {coeff}"
