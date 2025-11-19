# scf/scf_loop.py
# -*- coding: utf-8 -*-
"""
scf_loop
========

Este módulo contiene las rutinas principales del ciclo de autoconsistencia (SCF)
para el software **DFToy** en 1D. Incluye la generación de la densidad inicial,
la construcción del potencial efectivo y la iteración SCF hasta convergencia.

Contiene tres componentes principales:

- `InitialDensityParams`: clase que almacena los parámetros de la densidad inicial.
- `InitialDensityBuilder`: construye la densidad inicial a partir de un `Grid1D`.
- `SCFRunner`: ejecuta el ciclo SCF, calcula el estado fundamental y la energía total.

Ejemplo
-------
>>> from scf.scf_loop import InitialDensityParams, InitialDensityBuilder, SCFRunner
>>> from config import Config
>>> from grid import Grid1D
>>> from potential import ExternalPotentialParams, Potentials
>>> grid = Grid1D(Config.N, Config.x_min, Config.x_max)
>>> init_params = InitialDensityParams(guess_type=2, Np=2, sigma=1.0, offset=1.0)
>>> rho0 = InitialDensityBuilder.build(grid, init_params)
>>> ext_params = ExternalPotentialParams(potential_type=3, omega_ext=1.0)
>>> runner = SCFRunner(Config, grid, ext_params, m=1.0, k=1.0, rho0=2.0)
>>> result = runner.run(rho0)
>>> print(result['E'], result['rho'])
"""

import numpy as np
from dataclasses import dataclass
from pathlib import Path

from .config import Config
from .grid import Grid1D
from .potential import ExternalPotentialParams, Potentials
from .hamiltonian import Hamiltonian1D
from .solver import Solver


@dataclass
class InitialDensityParams:
    """
    Parámetros para construir la densidad inicial del ciclo SCF.

    Atributos
    ---------
    guess_type : int
        Tipo de densidad inicial:
        1 → gaussiana centrada en 0,
        2 → suma de dos gaussianas centradas en ±offset,
        otro → densidad tipo Lorentziana.
    Np : int
        Número de partículas o normalización de la densidad.
    sigma : float
        Ancho de las gaussianas iniciales.
    offset : float
        Desplazamiento de las gaussianas para `guess_type=2`.
    """
    guess_type: int
    Np: int
    sigma: float
    offset: float


class InitialDensityBuilder:
    """
    Construye la densidad inicial sobre la malla 1D.

    Métodos
    -------
    build(grid, init)
        Genera la densidad inicial normalizada a partir de un Grid1D y parámetros.
    """

    @staticmethod
    def build(grid: Grid1D, init: InitialDensityParams) -> np.ndarray:
        """
        Genera la densidad inicial :math:`ρ(x)` sobre la malla 1D.

        Parámetros
        ----------
        grid : Grid1D
            Malla espacial donde se define la densidad.
        init : InitialDensityParams
            Parámetros de la densidad inicial.

        Devuelve
        --------
        numpy.ndarray
            Arreglo de densidad inicial normalizada a `init.Np`.

        Notas
        -----
        Dependiendo de `init.guess_type`, se usa:
        - Gaussiana simple,
        - Suma de dos gaussianas desplazadas,
        - Densidad tipo Lorentziana.
        """
        x = grid.x
        if init.guess_type == 1:
            rho = np.exp(-(x * x) / (2.0 * init.sigma * init.sigma))
        elif init.guess_type == 2:
            rho = (np.exp(-((x - init.offset) ** 2) / (2.0 * init.sigma ** 2)) +
                   np.exp(-((x + init.offset) ** 2) / (2.0 * init.sigma ** 2)))
        else:
            rho = 1.0 / (1.0 + (x * x) / (init.sigma * init.sigma))

        integral = np.sum(rho) * grid.dx
        if integral > 0:
            rho = rho * (init.Np / integral)
        return rho


class SCFRunner:
    """
    Ejecuta el ciclo de autoconsistencia SCF.

    Atributos
    ---------
    cfg : Config
        Configuración global del cálculo.
    grid : Grid1D
        Malla espacial 1D.
    ext : ExternalPotentialParams
        Parámetros del potencial externo.
    m : float
        Masa de la partícula.
    k : float
        Constante del acoplamiento (interacción efectiva).
    rho0 : float
        Densidad inicial de referencia.

    Métodos
    -------
    run(rho_init, cancel_cb=None, q=None)
        Ejecuta el ciclo SCF hasta convergencia o cancelación.
    """

    def __init__(self, cfg: Config, grid: Grid1D, ext: ExternalPotentialParams, m: float, k: float, rho0: float):
        self.cfg = cfg
        self.grid = grid
        self.ext = ext
        self.m = m
        self.k = k
        self.rho0 = rho0

    def run(self, rho_init: np.ndarray, cancel_cb=None, q=None):
        """
        Ejecuta la iteración SCF hasta convergencia.

        Parámetros
        ----------
        rho_init : numpy.ndarray
            Densidad inicial sobre la malla.
        cancel_cb : callable, opcional
            Función que devuelve True para cancelar la ejecución.
        q : queue.Queue, opcional
            Cola para enviar mensajes de progreso/log.

        Devuelve
        --------
        dict
            Diccionario con:
            - "mu": energía del orbital más bajo,
            - "E": energía total del sistema,
            - "rho": densidad final convergida.

        Notas
        -----
        Calcula el potencial efectivo usando `Potentials.v_eff`,
        diagonaliza el Hamiltoniano tridiagonal con `Solver.lowest_eigenpair_tridiagonal`,
        aplica mezcla lineal de densidades y verifica convergencia según `cfg.tolerance`.
        """
        x = self.grid.x
        dx = self.grid.dx
        rho = rho_init.copy()
        conv = open("convergencia.dat", "w")

        omega0, omega = Potentials.compute_omegas(self.cfg.hbar, self.k, self.m, self.rho0)
        tcoeff = Config.kinetic_coeff(self.cfg.hbar, self.m, dx)  # negativo
        off = tcoeff * np.ones(self.grid.N - 1)

        for it in range(self.cfg.max_iter):
            if cancel_cb and cancel_cb():
                conv.close()
                if q:
                    q.put({"type": "log", "text": f"Cancelado en la iteración {it}."})
                return {"canceled": True}

            vext = Potentials.v_ext(x, self.ext, self.m)
            veff = Potentials.v_eff(x, vext, self.cfg.hbar, self.m, omega0, omega)

            diag = veff + (-2.0 * tcoeff)
            mu, phi0 = Solver.lowest_eigenpair_tridiagonal(diag, off)
            phi0 = Solver.normalize_phi(phi0, dx)
            rho_new = 2.0 * (phi0 * phi0)

            # Mezcla
            rho = (1.0 - self.cfg.mix_alpha) * rho + self.cfg.mix_alpha * rho_new

            # Error
            err = float(np.sqrt(np.sum((rho_new - rho)**2) * dx))
            conv.write(f"{it} {err}\n")
            if q and (it % 5 == 0):
                q.put({"type": "log", "text": f"Iter {it}: err={err:.2e}"})
            if err < self.cfg.tolerance:
                break

        E = 2.0 * mu - self.cfg.hbar * (omega0 + omega) * omega0 / (2.0 * omega)
        conv.close()
        return {"mu": mu, "E": E, "rho": rho}
