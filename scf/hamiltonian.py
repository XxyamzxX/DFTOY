# -*- coding: utf-8 -*-
import numpy as np
from dataclasses import dataclass

@dataclass
class Hamiltonian1D:
    diag: np.ndarray
    off: np.ndarray

    @staticmethod
    def build(diag: np.ndarray, off: np.ndarray):
        return Hamiltonian1D(diag=diag, off=off)

    def to_dense(self) -> np.ndarray:
        N = self.diag.size
        H = np.diag(self.diag) + np.diag(self.off, k=1) + np.diag(self.off, k=-1)
        return H
