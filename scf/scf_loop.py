# -*- coding: utf-8 -*-
import numpy as np
import time
from pathlib import Path
from dataclasses import dataclass

from .config import Config
from .grid import Grid1D
from .potential import ExternalPotentialParams, Potentials
from .hamiltonian import Hamiltonian1D
from .solver import Solver
from .io_utils import IOUtils

@dataclass
class InitialDensityParams:
    guess_type: int
    Np: int
    sigma: float
    offset: float

class InitialDensityBuilder:
    @staticmethod
    def build(grid: Grid1D, init: InitialDensityParams) -> np.ndarray:
        x = grid.x
        rho = np.zeros_like(x)
        if init.guess_type == 1:
            rho = np.exp(-(x*x) / (2.0 * init.sigma * init.sigma))
        elif init.guess_type == 2:
            rho = np.exp(-((x - init.offset)*(x - init.offset)) / (2.0 * init.sigma * init.sigma)) \
                + np.exp(-((x + init.offset)*(x + init.offset)) / (2.0 * init.sigma * init.sigma))
        elif init.guess_type == 3:
            rho = 1.0 / (1.0 + (x*x) / (init.sigma * init.sigma))
        integral = np.sum(rho) * grid.dx
        if integral <= 0:
            integral = 1e-12
        rho *= float(init.Np) / integral
        return rho

class SCFRunner:
    def __init__(self, cfg: Config, grid: Grid1D, ext_params: ExternalPotentialParams, m: float, k: float, rho0: float):
        self.cfg = cfg
        self.grid = grid
        self.ext_params = ext_params
        self.m = m
        self.k = k
        self.rho0 = rho0
        self.hbar = cfg.hbar

    def run(self, rho_init: np.ndarray):
        rho = rho_init.copy()
        rho_old = rho.copy()
        total_time = 0.0
        mu = 0.0
        omega0 = 0.0
        omega = 0.0

        conv_file = Path("convergencia.dat").open("w")
        for it in range(1, self.cfg.max_iter + 1):
            start = time.perf_counter()

            omega0, omega = Potentials.compute_omegas(self.hbar, self.k, self.m, self.rho0)
            vext = Potentials.v_ext(self.grid.x, self.ext_params, self.m)
            veff = Potentials.v_eff(self.grid.x, vext, self.hbar, self.m, omega0, omega)

            kin_coeff = self.cfg.kinetic_coeff(self.hbar, self.m, self.grid.dx)
            diag = (-2.0 * kin_coeff) + veff
            off = np.full(self.grid.N - 1, kin_coeff, dtype=float)

            mu, phi = Solver.lowest_eigenpair_tridiagonal(diag, off)
            phi = Solver.normalize_phi(phi, self.grid.dx)

            rho_new = 2.0 * (phi * phi)
            rho = self.cfg.mix_alpha * rho_new + (1.0 - self.cfg.mix_alpha) * rho_old

            elapsed = time.perf_counter() - start
            print(f"Iteración {it} tardó {elapsed:.6f} segundos.")
            total_time += elapsed

            diff = np.linalg.norm(rho - rho_old)
            conv_file.write(f"{it} {diff}\n")

            if diff < self.cfg.tolerance:
                print(f"Convergencia alcanzada en {it} iteraciones.")
                break

            rho_old = rho.copy()
            self.rho0 = float(np.max(rho))

        conv_file.close()

        E = 2.0 * mu - self.hbar * (omega0 + omega) * omega0 / (2.0 * omega)
        return mu, E, total_time, rho
