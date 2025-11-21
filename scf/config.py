# -*- coding: utf-8 -*-
"""
config
======

Define los parámetros globales y constantes físicas utilizadas por el paquete **DFToy**.

Este módulo centraliza todas las configuraciones numéricas y físicas del cálculo
—tamaño de malla, límites espaciales, parámetros de mezcla, valores iniciales, etc.—
para asegurar consistencia en todo el programa.

Contiene además dos métodos utilitarios que permiten calcular:
- el espaciamiento espacial `dx` de la malla, y  
- el coeficiente cinético para el operador laplaciano discreto.

Ejemplo
-------
>>> from DFToy.config import Config
>>> dx_value = Config.dx(Config.N, Config.x_min, Config.x_max)
>>> coeff = Config.kinetic_coeff(Config.hbar, m=1.0, dx=dx_value)
>>> print(dx_value, coeff)
"""

import math

class Config:
    """
    Clase contenedora de parámetros numéricos y físicos usados en DFToy.

    Esta clase funciona como un espacio de nombres donde se almacenan valores
    constantes que controlan la simulación: tamaño de la malla, límites
    espaciales, tolerancias, parámetros físicos como :math:`\\hbar`,
    parámetros de inicialización de densidad, y el coeficiente de mezcla
    utilizado en el ciclo de autoconsistencia (SCF).

    Atributos
    ---------
    N : int
        Número de puntos en la malla espacial 1D.
    x_min : float
        Límite izquierdo del dominio.
    x_max : float
        Límite derecho del dominio.
    tolerance : float
        Criterio de convergencia para el algoritmo SCF.
    max_iter : int
        Número máximo de iteraciones permitidas en el ciclo SCF.
    hbar : float
        Constante reducida de Planck usada en el cálculo del operador cinético.
    Np : int
        Número de densidades gaussianas iniciales usadas en el arranque.
    sigma : float
        Ancho de las gaussianas iniciales.
    offset : float
        Separación entre los centros de las gaussianas iniciales.
    mix_alpha : float
        Parámetro de mezcla lineal en el ciclo autoconsistente.

    Notas
    -----
    No es necesario instanciar esta clase; sus atributos se usan directamente
    mediante `Config.Atributo`.
    """

    # Parámetros numéricos
    N: int = 200
    x_min: float = -5.0
    x_max: float = 5.0
    tolerance: float = 1e-6
    max_iter: int = 100

    # Parámetros físicos
    hbar: float = 1.0

    # Parámetros de inicialización de densidad
    Np: int = 2
    sigma: float = 1.0
    offset: float = 1.0

    # Parámetro de mezcla en el SCF
    mix_alpha: float = 0.3

    @staticmethod
    def dx(N: int, x_min: float, x_max: float) -> float:
        """
        Calcula el espaciamiento espacial :math:`\\Delta x` de la malla 1D.

        Parámetros
        ----------
        N : int
            Número de puntos de la malla.
        x_min : float
            Límite izquierdo del dominio.
        x_max : float
            Límite derecho del dominio.

        Devuelve
        --------
        float
            Valor del espaciamiento :math:`\\Delta x`.

        Ejemplo
        -------
        >>> dx = Config.dx(200, -5.0, 5.0)
        >>> print(dx)
        0.05025125628140704
        """
        return (x_max - x_min) / (N - 1)

    @staticmethod
    def kinetic_coeff(hbar: float, m: float, dx: float) -> float:
        """
        Calcula el coeficiente cinético del operador laplaciano discreto.

        El coeficiente corresponde al término:

        .. math::

            - \\frac{\\hbar^2}{2 m \\, (\\Delta x)^2}

        Parámetros
        ----------
        hbar : float
            Constante reducida de Planck.
        m : float
            Masa efectiva del sistema.
        dx : float
            Espaciamiento espacial de la malla.

        Devuelve
        --------
        float
            Valor del coeficiente cinético.

        Ejemplo
        -------
        >>> dx = Config.dx(200, -5, 5)
        >>> coeff = Config.kinetic_coeff(1.0, m=1.0, dx=dx)
        >>> print(coeff)
        -198.0  # (aprox)
        """
        return -hbar * hbar / (2.0 * m * dx * dx)
