# -*- coding: utf-8 -*-
import math

class Config:
    # Parámetros numéricos
    N: int = 200
    x_min: float = -5.0
    x_max: float = 5.0
    tolerance: float = 1e-6
    max_iter: int = 100

    # Físicas
    hbar: float = 1.0

    # Inicialización de densidad
    Np: int = 2
    sigma: float = 1.0
    offset: float = 1.0

    # Mezcla
    mix_alpha: float = 0.3

    @staticmethod
    def dx(N: int, x_min: float, x_max: float) -> float:
        return (x_max - x_min) / (N - 1)

    @staticmethod
    def kinetic_coeff(hbar: float, m: float, dx: float) -> float:
        return -hbar * hbar / (2.0 * m * dx * dx)
