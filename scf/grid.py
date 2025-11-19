# -*- coding: utf-8 -*-
"""
grid
====

Define estructuras de malla espacial utilizadas en DFToy.  
Actualmente incluye la clase :class:`Grid1D`, una malla uniforme en una dimensión,
usada como base para representar funciones dependientes del espacio, operadores
discretizados y densidades electrónicas.

El módulo permite crear rápidamente una malla 1D especificando el número de nodos
y los límites del dominio.

Ejemplo
-------
>>> from DFToy.grid import Grid1D
>>> grid = Grid1D(N=200, x_min=-5.0, x_max=5.0)
>>> print(grid.dx)
0.05025125628140704
>>> print(grid.x[:5])
[-5.         -4.94974874 -4.89949749 -4.84924623 -4.79899497]
"""

import numpy as np
from dataclasses import dataclass


@dataclass
class Grid1D:
    """
    Malla espacial 1D uniforme para cálculos numéricos en DFToy.

    Esta clase genera una malla de puntos igualmente espaciados en el intervalo
    :math:`[x_{min}, x_{max}]`. La malla se usa para discretizar funciones,
    definir operadores diferenciales y evaluar densidades electrónicas.

    Parámetros
    ----------
    N : int
        Número total de puntos de la malla.
    x_min : float
        Extremo izquierdo del dominio espacial.
    x_max : float
        Extremo derecho del dominio espacial.

    Atributos
    ---------
    x : numpy.ndarray of shape (N,)
        Arreglo con las posiciones espaciales de cada nodo.
    dx : float
        Tamaño del espaciamiento :math:`\\Delta x` entre nodos consecutivos.

    Notas
    -----
    Esta clase utiliza ``dataclasses`` y calcula automáticamente los valores
    de ``x`` y ``dx`` durante la inicialización vía ``__post_init__``.
    """

    N: int
    x_min: float
    x_max: float

    def __post_init__(self):
        """
        Inicializa la malla espacial y calcula `dx`.

        Este método genera el arreglo de posiciones ``x`` usando
        :func:`numpy.linspace` y asigna el valor del espaciamiento ``dx``.

        Ejemplo
        -------
        >>> grid = Grid1D(N=5, x_min=0.0, x_max=1.0)
        >>> grid.x
        array([0.  , 0.25, 0.5 , 0.75, 1.  ])
        >>> grid.dx
        0.25
        """
        self.x = np.linspace(self.x_min, self.x_max, self.N)
        self.dx = (self.x_max - self.x_min) / (self.N - 1)
