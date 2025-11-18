# scf/scf_loop.py
# -*- coding: utf-8 -*-
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
    guess_type: int
    Np: int
    sigma: float
    offset: float

class InitialDensityBuilder:
    @staticmethod
    def build(grid: Grid1D, init: InitialDensityParams) -> np.ndarray:
        x = grid.x
        if init.guess_type == 1:
            rho = np.exp(-(x*x) / (2.0 * init.sigma * init.sigma))
        elif init.guess_type == 2:
            rho = (np.exp(-((x - init.offset)**2) / (2.0*init.sigma**2)) +
                   np.exp(-((x + init.offset)**2) / (2.0*init.sigma**2)))
        else:
            rho = 1.0 / (1.0 + (x*x) / (init.sigma * init.sigma))
        # Normalizar a Np
        integral = np.sum(rho) * grid.dx
        if integral > 0:
            rho = rho * (init.Np / integral)
        return rho

class SCFRunner:
    def __init__(self, cfg: Config, grid: Grid1D, ext: ExternalPotentialParams, m: float, k: float, rho0: float):
        self.cfg = cfg
        self.grid = grid
        self.ext = ext
        self.m = m
        self.k = k
        self.rho0 = rho0

    def run(self, rho_init: np.ndarray, cancel_cb=None, q=None):
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
                if q: q.put({"type":"log", "text": f"Cancelado en la iteración {it}."})
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
                q.put({"type":"log", "text": f"Iter {it}: err={err:.2e}"})
            if err < self.cfg.tolerance:
                break

        E = 2.0 * mu - self.cfg.hbar * (omega0 + omega) * omega0 / (2.0 * omega)
        conv.close()
        return {"mu": mu, "E": E, "rho": rho}
