# scf/solver.py
# -*- coding: utf-8 -*-
"""
solver
======

Proporciona métodos para resolver problemas de valores y vectores propios
de matrices tridiagonales 1D, típicamente usadas para el Hamiltoniano
discreto en DFToy. Incluye normalización de funciones de onda.

Contiene la clase:
- `Solver`: métodos estáticos para cálculo de eigenpares y normalización.
  
Ejemplo
-------
>>> import numpy as np
>>> from scf.solver import Solver
>>> diag = np.array([1.0, 2.0, 3.0])
>>> off = np.array([-1.0, -1.0])
>>> mu, phi0 = Solver.lowest_eigenpair_tridiagonal(diag, off)
>>> print(mu)
"""

import numpy as np
import math


class Solver:
    """
    Métodos estáticos para resolver eigenpares de matrices tridiagonales
    y normalizar funciones de onda.
    """

    @staticmethod
    def lowest_eigenpair_tridiagonal(diag: np.ndarray, off: np.ndarray):
        """
        Calcula el valor propio más bajo y su vector propio asociado
        de una matriz tridiagonal.

        Parámetros
        ----------
        diag : np.ndarray
            Elementos de la diagonal principal.
        off : np.ndarray
            Elementos de las diagonales off-diagonal (k=±1).

        Devuelve
        -------
        tuple
            (valor propio más bajo, vector propio correspondiente).

        Ejemplo
        -------
        >>> diag = np.array([2.0, 3.0, 4.0])
        >>> off = np.array([-1.0, -1.0])
        >>> mu, phi0 = Solver.lowest_eigenpair_tridiagonal(diag, off)
        """
        H = np.diag(diag) + np.diag(off, k=1) + np.diag(off, k=-1)
        w, V = np.linalg.eigh(H)
        return w[0], V[:, 0]

    @staticmethod
    def eigenpairs_tridiagonal(diag: np.ndarray, off: np.ndarray, k: int):
        """
        Calcula los k valores y vectores propios más bajos de una matriz tridiagonal.

        Parámetros
        ----------
        diag : np.ndarray
            Elementos de la diagonal principal.
        off : np.ndarray
            Elementos de las diagonales off-diagonal (k=±1).
        k : int
            Número de eigenpares a devolver.

        Devuelve
        -------
        tuple
            (valores propios, vectores propios) de los k estados más bajos.

        Ejemplo
        -------
        >>> diag = np.array([2.0, 3.0, 4.0])
        >>> off = np.array([-1.0, -1.0])
        >>> w, V = Solver.eigenpairs_tridiagonal(diag, off, 2)
        """
        H = np.diag(diag) + np.diag(off, k=1) + np.diag(off, k=-1)
        w, V = np.linalg.eigh(H)
        k = min(k, len(w))
        return w[:k], V[:, :k]

    @staticmethod
    def normalize_phi(phi: np.ndarray, dx: float) -> np.ndarray:
        """
        Normaliza una función de onda discreta según la integral discreta
        sobre la malla (sumatoria multiplicada por dx).

        Parámetros
        ----------
        phi : np.ndarray
            Vector de la función de onda.
        dx : float
            Espaciamiento de la malla.

        Devuelve
        -------
        np.ndarray
            Función de onda normalizada.

        Lanza
        -----
        RuntimeError
            Si la norma calculada no es positiva.

        Ejemplo
        -------
        >>> phi = np.array([1.0, 2.0, 1.0])
        >>> phi_norm = Solver.normalize_phi(phi, dx=0.1)
        """
        norm_phi = math.sqrt(float(np.sum(phi * phi) * dx))
        if norm_phi <= 0:
            raise RuntimeError("Norma no positiva.")
        return phi / norm_phi

