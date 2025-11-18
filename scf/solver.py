# scf/solver.py
# -*- coding: utf-8 -*-
import numpy as np
import math

class Solver:
    @staticmethod
    def lowest_eigenpair_tridiagonal(diag: np.ndarray, off: np.ndarray):
        H = np.diag(diag) + np.diag(off, k=1) + np.diag(off, k=-1)
        w, V = np.linalg.eigh(H)
        return w[0], V[:, 0]

    @staticmethod
    def eigenpairs_tridiagonal(diag: np.ndarray, off: np.ndarray, k: int):
        H = np.diag(diag) + np.diag(off, k=1) + np.diag(off, k=-1)
        w, V = np.linalg.eigh(H)
        k = min(k, len(w))
        return w[:k], V[:, :k]

    @staticmethod
    def normalize_phi(phi: np.ndarray, dx: float) -> np.ndarray:
        norm_phi = math.sqrt(float(np.sum(phi * phi) * dx))
        if norm_phi <= 0:
            raise RuntimeError("Norma no positiva.")
        return phi / norm_phi
