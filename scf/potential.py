# -*- coding: utf-8 -*-
import numpy as np
import math
from dataclasses import dataclass

@dataclass
class ExternalPotentialParams:
    potential_type: int
    alpha: float = 1.0
    beta: float = 1.0
    a_dw: float = 1.0
    omega_ext: float = 1.0
    V0: float = 1.0
    L: float = 2.0

class Potentials:
    @staticmethod
    def v_ext(x: np.ndarray, p: ExternalPotentialParams, m: float) -> np.ndarray:
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
        a_term = m * (2.0 * ((omega * omega0) / (omega0 + omega))**2 - 0.5 * omega0 * omega0)
        b_term = hbar * ((omega0 + omega)**2) / (4.0 * omega)
        c_term = hbar * (omega * omega0) / (omega0 + omega)
        return v_ext + a_term * (x*x) + b_term - c_term
