# -*- coding: utf-8 -*-
import time
import numpy as np
from pathlib import Path

from scf.config import Config
from scf.grid import Grid1D
from scf.potential import ExternalPotentialParams
from scf.scf_loop import InitialDensityParams, InitialDensityBuilder, SCFRunner
from scf.io_utils import IOUtils

def run_scf_once(params: dict, q):
    """
    Ejecuta el SCF con los parámetros actuales y comunica progreso/resultado
    a la GUI usando la cola q (mensajes 'log', 'done' o 'error').
    """
    try:
        cfg = Config()
        # Inyectar tolerancia, etc., si quisieras desde GUI
        grid = Grid1D(cfg.N, cfg.x_min, cfg.x_max)

        init = InitialDensityParams(
            guess_type=int(params.get("guess_type", 1)),
            Np=cfg.Np,
            sigma=float(params.get("sigma", 1.0)),
            offset=float(params.get("offset", 1.0)),
        )
        rho_init = InitialDensityBuilder.build(grid, init)

        ext = ExternalPotentialParams(
            potential_type=int(params.get("potential_type", 3)),
            alpha=float(params.get("alpha", 1.0)),
            beta=float(params.get("beta", 1.0)),
            a_dw=float(params.get("a_dw", 1.0)),
            omega_ext=float(params.get("omega_ext", 1.0)),
            V0=float(params.get("V0", 1.0)),
            L=float(params.get("L", 2.0)),
        )

        k = float(params.get("k", 1.0))
        m = float(params.get("m", 1.0))
        rho0 = float(params.get("rho0", 1.0))

        q.put({"kind": "log", "text": "Preparando datos y malla..."})
        t0 = time.perf_counter()
        runner = SCFRunner(cfg, grid, ext, m=m, k=k, rho0=rho0)

        mu, E, total_time, rho = runner.run(rho_init)

        # Guardar y crear scripts/plots como en el CLI
        IOUtils.save_density(grid.x, rho, "rho_vs_x.dat")
        IOUtils.write_gnuplot_density_script(E, int(params.get("potential_type", 3)), "plot.gp")
        IOUtils.write_gnuplot_convergence_script("convergencia.gp")
        # Los EPS se generan después en el hilo UI llamando IOUtils.run_gnuplot

        t1 = time.perf_counter()
        q.put({"kind": "done", "mu": mu, "E": E, "time": t1 - t0})
    except Exception as e:
        q.put({"kind": "error", "text": f"{type(e).__name__}: {e}"})
