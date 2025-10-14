# -*- coding: utf-8 -*-
import numpy as np
from dataclasses import dataclass

@dataclass
class Grid1D:
    N: int
    x_min: float
    x_max: float

    def __post_init__(self):
        self.x = np.linspace(self.x_min, self.x_max, self.N)
        self.dx = (self.x_max - self.x_min) / (self.N - 1)