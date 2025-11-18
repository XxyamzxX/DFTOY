# gui_app/bridge.py
# -*- coding: utf-8 -*-
import time
from pathlib import Path
import numpy as np

from scf.config import Config
from scf.grid import Grid1D
from scf.potential import ExternalPotentialParams
from scf.scf_loop import InitialDensityParams, InitialDensityBuilder, SCFRunner
from scf.io_utils import IOUtils

def run_scf_once(params: dict, q, cancel_evt=None):
    try:
        cfg = Config()
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

        q.put({"type":"log", "text":"Preparando datos y malla..."})
        cancel_cb = (lambda: cancel_evt.is_set()) if cancel_evt is not None else (lambda: False)

        t0 = time.perf_counter()
        runner = SCFRunner(cfg, grid, ext, m=m, k=k, rho0=rho0)
        status = runner.run(rho_init, cancel_cb=cancel_cb, q=q)
        if status.get("canceled"):
            q.put({"type":"canceled"})
            return
        mu = status["mu"]; E = status["E"]; rho = status["rho"]

        IOUtils.save_density(grid.x, rho, "rho_vs_x.dat")
        IOUtils.write_gnuplot_density_script(E, int(params.get("potential_type", 3)), "plot.gp")
        IOUtils.write_gnuplot_convergence_script("convergencia.gp")

        t1 = time.perf_counter()
        q.put({"type":"done", "mu": mu, "E": E, "time": t1 - t0})
    except Exception as e:
        q.put({"type":"error", "text": f"{type(e).__name__}: {e}"})
