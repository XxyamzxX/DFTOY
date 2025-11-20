# tests/test_solver.py
# -*- coding: utf-8 -*-
"""
Unit tests para solver.py

POR QUÉ: Solver calcula valores y vectores propios del Hamiltoniano.
Errores aquí causan energías fundamentales incorrectas, convergencia
falsa y funciones de onda no normalizadas.
"""

import pytest
import numpy as np
from scf.solver import Solver

def test_lowest_eigenpair_simple():
    """
    Verifica cálculo de eigenpar más bajo en matriz simple.
    
    POR QUÉ: Para un Hamiltoniano simple H=diag([1,2,3]) con off=0,
    el menor eigenvalor debe ser 1. Detecta errores en np.linalg.eigh
    o indexación w[0].
    """
    diag = np.array([3.0, 2.0, 1.0])
    off = np.zeros(2)
    mu, phi = Solver.lowest_eigenpair_tridiagonal(diag, off)
    assert np.isclose(mu, 1.0, atol=1e-9), f"Eigenvalor incorrecto: {mu}"
    assert phi.shape == (3,), "Eigenvector debe tener shape (N,)"

def test_normalize_phi_integral():
    """
    Verifica que normalize_phi deje ∫|φ|²dx = 1.
    
    POR QUÉ: La normalización es requisito físico; sin ella, densidades
    y probabilidades son incorrectas. Suma(φ²·dx) debe ser ~1 tras normalizar.
    """
    phi = np.array([1.0, 2.0, 1.0])
    dx = 0.1
    phi_norm = Solver.normalize_phi(phi, dx)
    integral = np.sum(phi_norm**2) * dx
    assert np.isclose(integral, 1.0, rtol=1e-6), f"Norma incorrecta: {integral}"

def test_normalize_phi_raises_on_zero():
    """
    Verifica que normalize_phi lance error con φ=0.
    
    POR QUÉ: Un vector nulo no tiene norma definida; el programa
    debe fallar explícitamente en lugar de dividir por cero.
    """
    phi = np.zeros(5)
    dx = 0.1
    with pytest.raises(RuntimeError, match="Norma no positiva"):
        Solver.normalize_phi(phi, dx)

def test_eigenpairs_tridiagonal_count():
    """
    Verifica que eigenpairs_tridiagonal retorne k eigenpares.
    
    POR QUÉ: Para graficar estados excitados se necesitan múltiples eigenpares;
    si k > tamaño de la matriz, debe truncarse sin crash.
    """
    diag = np.array([4.0, 3.0, 2.0, 1.0])
    off = np.array([-0.5, -0.5, -0.5])
    w, V = Solver.eigenpairs_tridiagonal(diag, off, k=2)
    assert w.shape == (2,), f"Debe devolver 2 eigenvalores, obtuvo {w.shape}"
    assert V.shape == (4, 2), f"Debe devolver matriz (4,2), obtuvo {V.shape}"

def test_eigenpairs_sorted():
    """
    Verifica que eigenvalores estén ordenados ascendentemente.
    
    POR QUÉ: np.linalg.eigh garantiza orden ascendente; si no, el "estado
    fundamental" puede ser un estado excitado, rompiendo todo el SCF.
    """
    diag = np.array([5.0, 1.0, 3.0])
    off = np.zeros(2)
    w, _ = Solver.eigenpairs_tridiagonal(diag, off, k=3)
    assert np.all(w[:-1] <= w[1:]), "Eigenvalores deben estar ordenados"
