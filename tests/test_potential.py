# tests/test_potential.py
# -*- coding: utf-8 -*-
"""
Unit tests para potential.py

POR QUÉ: Los potenciales definen la física del problema. Errores en v_ext
o v_eff generan estados fundamentales incorrectos, energías sin sentido
y fallas de convergencia en el SCF.
"""

import pytest
import numpy as np
from scf.potential import ExternalPotentialParams, Potentials

def test_v_ext_linear():
    """
    Verifica potencial lineal V(x) = α|x|.
    
    POR QUÉ: Debe ser simétrico, cero en x=0 y crecer linealmente.
    Detecta errores en condicional potential_type==1.
    """
    x = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    p = ExternalPotentialParams(potential_type=1, alpha=2.0)
    v = Potentials.v_ext(x, p, m=1.0)
    expected = np.array([4.0, 2.0, 0.0, 2.0, 4.0])
    assert np.allclose(v, expected), f"v_ext lineal incorrecto: {v}"

def test_v_ext_double_well():
    """
    Verifica potencial doble pozo V(x) = β(x²-a²)².
    
    POR QUÉ: Debe tener mínimos en x=±a y máximo en x=0.
    Detecta errores de signo o paréntesis en la fórmula cuártica.
    """
    x = np.array([-1.0, 0.0, 1.0])
    p = ExternalPotentialParams(potential_type=2, beta=1.0, a_dw=1.0)
    v = Potentials.v_ext(x, p, m=1.0)
    assert v[1] > v[0], "Centro debe ser máximo local"
    assert v[1] > v[2], "Centro debe ser máximo local"
    assert np.isclose(v[0], 0.0, atol=1e-9), "Mínimo en x=-a debe ser ~0"
    assert np.isclose(v[2], 0.0, atol=1e-9), "Mínimo en x=+a debe ser ~0"

def test_v_ext_harmonic():
    """
    Verifica potencial armónico V(x) = 0.5·m·ω²·x².
    
    POR QUÉ: Debe ser cuadrático con mínimo en x=0.
    Detecta errores en coeficiente 0.5 o uso de omega_ext.
    """
    x = np.linspace(-2, 2, 5)
    p = ExternalPotentialParams(potential_type=3, omega_ext=2.0)
    m = 1.5
    v = Potentials.v_ext(x, p, m)
    expected = 0.5 * m * (2.0**2) * (x**2)
    assert np.allclose(v, expected, atol=1e-9), "Armónico incorrecto"

def test_v_ext_square_well():
    """
    Verifica pozo cuadrado: V(x) = -V0 si |x|<L/2, sino 0.
    
    POR QUÉ: Debe ser discontinuo en bordes ±L/2.
    Detecta errores en np.where o condición de |x|.
    """
    x = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])
    p = ExternalPotentialParams(potential_type=4, V0=3.0, L=1.5)
    v = Potentials.v_ext(x, p, m=1.0)
    # |x|<L/2=0.75 → índices 1,2,3
    assert v[0] == 0.0, "Fuera del pozo debe ser 0"
    assert v[1] == -3.0, "Dentro debe ser -V0"
    assert v[2] == -3.0
    assert v[3] == -3.0
    assert v[4] == 0.0

def test_compute_omegas_positive():
    """
    Verifica que compute_omegas retorne valores positivos.
    
    POR QUÉ: omega0 y omega son frecuencias físicas; si son negativas
    o NaN, la fórmula de v_eff falla. Detecta divisiones por cero
    o raíces de números negativos.
    """
    omega0, omega = Potentials.compute_omegas(hbar=1.0, k=1.0, m=1.0, rho0=1.0)
    assert omega0 > 0, "omega0 debe ser positiva"
    assert omega > 0, "omega debe ser positiva"
    assert not np.isnan(omega0), "omega0 no debe ser NaN"
    assert not np.isnan(omega), "omega no debe ser NaN"

def test_v_eff_shape():
    """
    Verifica que v_eff devuelva array del mismo tamaño que x.
    
    POR QUÉ: v_eff combina v_ext con términos cuadráticos y constantes;
    broadcasting erróneo causa shape mismatch y crash del SCF.
    """
    x = np.linspace(-3, 3, 50)
    v_ext = np.ones_like(x)
    v = Potentials.v_eff(x, v_ext, hbar=1.0, m=1.0, omega0=0.5, omega=1.0)
    assert v.shape == x.shape, "v_eff debe tener mismo shape que x"

def test_v_eff_adds_terms():
    """
    Verifica que v_eff = v_ext + término_cuadrático + constantes.
    
    POR QUÉ: v_eff debe modificar v_ext con términos específicos;
    si falta alguno, el potencial efectivo es incorrecto y SCF diverge.
    """
    x = np.array([0.0, 1.0])
    v_ext = np.zeros_like(x)
    v = Potentials.v_eff(x, v_ext, hbar=1.0, m=1.0, omega0=1.0, omega=2.0)
    # En x=0, v[0] debe ser constantes b_term - c_term
    # En x≠0, debe incluir a_term·x²
    assert v[1] != v[0], "v_eff debe incluir término cuadrático a·x²"
