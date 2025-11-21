# -*- coding: utf-8 -*-
"""
potentials
==========

Define los potenciales externos y efectivos utilizados en DFToy, junto con
funciones auxiliares para cálculos relacionados, como la determinación de
frecuencias características omega0 y omega.

Contiene:
- Clase `ExternalPotentialParams` para parametrizar distintos tipos de potenciales.
- Clase `Potentials` con métodos estáticos para calcular potencial externo, efectivo y frecuencias.

Ejemplo
-------
>>> import numpy as np
>>> from potentials import ExternalPotentialParams, Potentials
>>> x = np.linspace(-5,5,100)
>>> params = ExternalPotentialParams(potential_type=1, alpha=2.0)
>>> v = Potentials.v_ext(x, params, m=1.0)
>>> omega0, omega = Potentials.compute_omegas(hbar=1.0, k=1.0, m=1.0, rho0=1.0)
>>> v_eff_array = Potentials.v_eff(x, v, hbar=1.0, m=1.0, omega0=omega0, omega=omega)
"""

import numpy as np
import math
from dataclasses import dataclass

@dataclass
class ExternalPotentialParams:
    """
    Parámetros que definen un potencial externo.

    Atributos
    ---------
    potential_type : int
        Tipo de potencial:
        1 → V(x) = α|x|
        2 → V(x) = β(x² - a²)²
        3 → V(x) = 0.5 m ω² x²
        4 → Pozo cuadrado
    alpha : float, opcional
        Parámetro α usado en el potencial lineal. Por defecto 1.0
    beta : float, opcional
        Parámetro β usado en el potencial de doble pozo. Por defecto 1.0
    a_dw : float, opcional
        Parámetro a del doble pozo. Por defecto 1.0
    omega_ext : float, opcional
        Frecuencia para el potencial armónico. Por defecto 1.0
    V0 : float, opcional
        Profundidad del pozo cuadrado. Por defecto 1.0
    L : float, opcional
        Ancho del pozo cuadrado. Por defecto 2.0
    """

    potential_type: int
    alpha: float = 1.0
    beta: float = 1.0
    a_dw: float = 1.0
    omega_ext: float = 1.0
    V0: float = 1.0
    L: float = 2.0

class Potentials:
    """
    Métodos estáticos para calcular potenciales y frecuencias características.

    Contiene funciones para:
    - v_ext: calcular el potencial externo según el tipo de potencial.
    - compute_omegas: determinar las frecuencias omega0 y omega.
    - v_eff: calcular el potencial efectivo a partir del externo y las frecuencias.
    """

    @staticmethod
    def v_ext(x: np.ndarray, p: ExternalPotentialParams, m: float) -> np.ndarray:
        """
        Calcula el potencial externo V(x) según el tipo especificado.

        Parámetros
        ----------
        x : np.ndarray
            Vector de posiciones espaciales.
        p : ExternalPotentialParams
            Parámetros del potencial.
        m : float
            Masa de la partícula.

        Devuelve
        -------
        np.ndarray
            Vector con el valor del potencial externo en cada posición x.

        Ejemplo
        -------
        >>> import numpy as np
        >>> params = ExternalPotentialParams(potential_type=1, alpha=2.0)
        >>> x = np.array([-1,0,1])
        >>> Potentials.v_ext(x, params, m=1.0)
        array([2., 0., 2.])
        """
        if p.potential_type == 1:
            return p.alpha * np.abs(x)
        if p.potential_type == 2:
            return p.beta * (x*x - p.a_dw*p.a_dw)**2
        if p.potential_type == 3:
            return 0.5 * m * (p.omega_ext**2) * (x*x)
        if p.potential_type == 4:
            return np.where(np.abs(x) < (p.L/2.0), -p.V0, 0.0)
        return np.zeros_like(x)

    @staticmethod
    def compute_omegas(hbar: float, k: float, m: float, rho0: float):
        """
        Calcula las frecuencias características omega0 y omega usadas en el potencial efectivo.

        Parámetros
        ----------
        hbar : float
            Constante reducida de Planck.
        k : float
            Constante de fuerza.
        m : float
            Masa de la partícula.
        rho0 : float
            Valor característico de la densidad inicial.

        Devuelve
        -------
        tuple
            omega0, omega : float
            Frecuencias características calculadas.

        Notas
        -----
        La fórmula proviene de la aproximación para el potencial efectivo 1D.

        Ejemplo
        -------
        >>> omega0, omega = Potentials.compute_omegas(hbar=1.0, k=1.0, m=1.0, rho0=1.0)
        """
        c = (8.0 * m) / (math.pi * hbar * (rho0**2))
        c4 = c**4
        km = k / m
        gamma = 4.0 * k / (m * c4) \
              + (2.0 * k / (3.0 * m))**3 \
              + (4.0/9.0) * k / (m * c4) * math.sqrt(3.0 * (27.0 + 4.0 * c4 * (km**2)))
        gamma13 = gamma ** (1.0/3.0)
        y_term = gamma13 + ((2.0 * k) / (3.0 * m))**2 / gamma13 + (2.0 * k) / (3.0 * m)
        R = math.sqrt(1.0 / (c * c) - 2.0 * k / m + y_term)
        D = math.sqrt(3.0 / (c * c) - R * R - 4.0 * k / m + (4.0 * k / (m * c) + 2.0 / (c**3)) / R)
        omega0 = (1.0 / c + R + D) / 2.0
        omega = math.sqrt(omega0 * omega0 + 2.0 * k / m)
        return omega0, omega

    @staticmethod
    def v_eff(x: np.ndarray, v_ext: np.ndarray, hbar: float, m: float, omega0: float, omega: float) -> np.ndarray:
        """
        Calcula el potencial efectivo V_eff(x) a partir del potencial externo y las frecuencias.

        Parámetros
        ----------
        x : np.ndarray
            Vector de posiciones espaciales.
        v_ext : np.ndarray
            Vector con el potencial externo.
        hbar : float
            Constante reducida de Planck.
        m : float
            Masa de la partícula.
        omega0 : float
            Frecuencia característica omega0.
        omega : float
            Frecuencia característica omega.

        Devuelve
        -------
        np.ndarray
            Vector con el valor del potencial efectivo en cada posición x.

        Ejemplo
        -------
        >>> x = np.array([-1,0,1])
        >>> v_ext = np.array([1,0,1])
        >>> Potentials.v_eff(x, v_ext, hbar=1.0, m=1.0, omega0=0.5, omega=1.0)
        array([...])
        """
        a_term = m * (2.0 * ((omega * omega0) / (omega0 + omega))**2 - 0.5 * omega0 * omega0)
        b_term = hbar * ((omega0 + omega)**2) / (4.0 * omega)
        c_term = hbar * (omega * omega0) / (omega0 + omega)
        return v_ext + a_term * (x*x) + b_term - c_term
