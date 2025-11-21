# -*- coding: utf-8 -*-
"""
hamiltonian
===========

Define estructuras para representar el Hamiltoniano 1D en forma tridiagonal,
como se usa en métodos de diferencias finitas o discretización estándar en DFToy.

La clase :class:`Hamiltonian1D` contiene las diagonales necesarias para construir
la matriz completa y permite convertirla a una representación densa.

Ejemplo
-------
>>> import numpy as np
>>> from DFToy.hamiltonian import Hamiltonian1D
>>>
>>> diag = np.array([2., 2., 2.])
>>> off = np.array([-1., -1.])
>>> H = Hamiltonian1D.build(diag, off)
>>> H.to_dense()
array([[ 2., -1.,  0.],
       [-1.,  2., -1.],
       [ 0., -1.,  2.]])
"""

import numpy as np
from dataclasses import dataclass

@dataclass
class Hamiltonian1D:
    """
    Representación tridiagonal de un Hamiltoniano 1D.

    Esta clase almacena:
    - ``diag``: la diagonal principal del Hamiltoniano.
    - ``off``: la diagonal secundaria (superior e inferior), asumida simétrica.

    El Hamiltoniano tridiagonal típico para problemas 1D discretizados
    tiene la forma:

    .. math::

        H = \\begin{pmatrix}
        d_0 & o_0 & 0   & \\cdots & 0 \\\\
        o_0 & d_1 & o_1 & \\cdots & 0 \\\\
        0   & o_1 & d_2 & \\cdots & 0 \\\\
        \\vdots & \\vdots & \\vdots & \\ddots & \\vdots \\\\
        0 & 0 & 0 & \\cdots & d_{N-1}
        \\end{pmatrix}

    Parámetros
    ----------
    diag : numpy.ndarray
        Diagonal principal del Hamiltoniano (tamaño N).
    off : numpy.ndarray
        Diagonal secundaria superior/inferior (tamaño N-1).

    Notas
    -----
    Esta representación es eficiente para construir rápidamente el operador cinético
    o potencial en esquemas de diferencias finitas.
    """

    diag: np.ndarray
    off: np.ndarray

    @staticmethod
    def build(diag: np.ndarray, off: np.ndarray):
        """
        Construye una instancia de :class:`Hamiltonian1D` a partir de
        sus diagonales.

        Parámetros
        ----------
        diag : numpy.ndarray
            Arreglo 1D con la diagonal principal.
        off : numpy.ndarray
            Arreglo 1D con los elementos de la diagonal secundaria.

        Devuelve
        --------
        Hamiltonian1D
            Instancia que contiene las diagonales especificadas.

        Ejemplo
        -------
        >>> diag = np.array([1., 1., 1.])
        >>> off = np.array([-0.5, -0.5])
        >>> H = Hamiltonian1D.build(diag, off)
        """
        return Hamiltonian1D(diag=diag, off=off)

    def to_dense(self) -> np.ndarray:
        """
        Construye la matriz densa correspondiente al Hamiltoniano tridiagonal.

        Devuelve
        --------
        numpy.ndarray of shape (N, N)
            Matriz densa del Hamiltoniano.

        Ejemplo
        -------
        >>> diag = np.array([2., 2., 2.])
        >>> off = np.array([-1., -1.])
        >>> H = Hamiltonian1D(diag, off)
        >>> H.to_dense()
        array([[ 2., -1.,  0.],
               [-1.,  2., -1.],
               [ 0., -1.,  2.]])
        """
        N = self.diag.size
        H = np.diag(self.diag) + np.diag(self.off, k=1) + np.diag(self.off, k=-1)
        return H
