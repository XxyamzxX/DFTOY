# tests/test_scf_loop.py
# -*- coding: utf-8 -*-
"""
Unit tests para scf_loop.py

POR QUÉ: SCFRunner ejecuta el ciclo autoconsistente; errores aquí causan
no-convergencia, energías erróneas y salidas con shape incorrecto.
"""

import pytest
import numpy as np
from scf.config import Config
from scf.grid import Grid1D
from scf.potential import ExternalPotentialParams
from scf.scf_loop import InitialDensityParams, InitialDensityBuilder, SCFRunner

def test_initial_density_normalization():
    """
    Verifica que la densidad inicial se normalice a Np.
    
    POR QUÉ: Si ρ no está normalizada, el número de partículas es incorrecto
    y la física del sistema es errónea desde el inicio.
    """
    grid = Grid1D(N=201, x_min=-5.0, x_max=5.0)
    init = InitialDensityParams(guess_type=1, Np=2, sigma=1.0, offset=1.0)
    rho = InitialDensityBuilder.build(grid, init)
    integral = np.sum(rho) * grid.dx
    assert np.isclose(integral, init.Np, rtol=1e-3), f"Normalización incorrecta: {integral} vs {init.Np}"

def test_initial_density_bimodal_has_two_peaks():
    """
    Verifica que guess_type=2 (bimodal) tenga dos picos.
    
    POR QUÉ: La densidad bimodal suma dos gaussianas centradas en ±offset;
    si falta un término o hay error de signo, solo aparece un pico.
    """
    grid = Grid1D(N=201, x_min=-5.0, x_max=5.0)
    init = InitialDensityParams(guess_type=2, Np=2, sigma=0.5, offset=1.5)
    rho = InitialDensityBuilder.build(grid, init)
    # Buscar índices cerca de ±offset
    idx_neg = np.argmin(np.abs(grid.x + init.offset))
    idx_pos = np.argmin(np.abs(grid.x - init.offset))
    # rho debe tener máximos locales cerca de esos puntos
    assert rho[idx_neg] > 0.5 * rho.max(), "Falta pico en x=-offset"
    assert rho[idx_pos] > 0.5 * rho.max(), "Falta pico en x=+offset"

def test_scf_runner_returns_dict():
    """
    Verifica que SCFRunner.run retorne dict con 'mu', 'E', 'rho'.
    
    POR QUÉ: El bridge y GUI esperan este formato; si cambia,
    se rompe la interfaz entre módulos.
    """
    cfg = Config()
    grid = Grid1D(cfg.N, cfg.x_min, cfg.x_max)
    ext = ExternalPotentialParams(potential_type=3, omega_ext=1.0)
    init = InitialDensityParams(guess_type=1, Np=cfg.Np, sigma=1.0, offset=1.0)
    rho0 = InitialDensityBuilder.build(grid, init)
    runner = SCFRunner(cfg, grid, ext, m=1.0, k=1.0, rho0=1.0)
    result = runner.run(rho0)
    assert isinstance(result, dict), "run debe retornar un diccionario"
    assert "mu" in result, "Falta 'mu' en resultado"
    assert "E" in result, "Falta 'E' en resultado"
    assert "rho" in result, "Falta 'rho' en resultado"

def test_scf_runner_rho_shape():
    """
    Verifica que la densidad final tenga shape correcto.
    
    POR QUÉ: Si rho cambia de tamaño durante SCF, causa crashes
    al comparar con rho_new o al guardar en archivos.
    """
    cfg = Config()
    grid = Grid1D(cfg.N, cfg.x_min, cfg.x_max)
    ext = ExternalPotentialParams(potential_type=3, omega_ext=1.0)
    init = InitialDensityParams(guess_type=1, Np=2, sigma=1.0, offset=1.0)
    rho_init = InitialDensityBuilder.build(grid, init)
    runner = SCFRunner(cfg, grid, ext, m=1.0, k=1.0, rho0=1.0)
    result = runner.run(rho_init)
    assert result["rho"].shape == (cfg.N,), f"Shape incorrecto: {result['rho'].shape}"

def test_scf_cancellation():
    """
    Verifica que SCF pueda cancelarse con cancel_cb.
    
    POR QUÉ: Si cancel_cb no se consulta correctamente, el botón
    "Cancelar" de la GUI no funciona y el programa se cuelga.
    """
    cfg = Config()
    cfg.max_iter = 100  # asegurar muchas iteraciones
    grid = Grid1D(cfg.N, cfg.x_min, cfg.x_max)
    ext = ExternalPotentialParams(potential_type=3, omega_ext=1.0)
    init = InitialDensityParams(guess_type=1, Np=2, sigma=1.0, offset=1.0)
    rho_init = InitialDensityBuilder.build(grid, init)
    runner = SCFRunner(cfg, grid, ext, m=1.0, k=1.0, rho0=1.0)
    
    cancel_cb = lambda: True  # cancelar inmediatamente
    result = runner.run(rho_init, cancel_cb=cancel_cb)
    assert result.get("canceled") is True, "SCF debe retornar canceled=True"

def test_scf_convergence_file_created():
    """
    Verifica que SCF cree archivo convergencia.dat.
    
    POR QUÉ: Este archivo registra el error por iteración; si falta,
    la GUI no puede graficar convergencia y se pierde diagnóstico.
    """
    import os
    cfg = Config()
    cfg.max_iter = 5  # pocas iteraciones
    grid = Grid1D(cfg.N, cfg.x_min, cfg.x_max)
    ext = ExternalPotentialParams(potential_type=3, omega_ext=1.0)
    init = InitialDensityParams(guess_type=1, Np=2, sigma=1.0, offset=1.0)
    rho_init = InitialDensityBuilder.build(grid, init)
    runner = SCFRunner(cfg, grid, ext, m=1.0, k=1.0, rho0=1.0)
    
    # Limpiar archivo previo
    if os.path.exists("convergencia.dat"):
        os.remove("convergencia.dat")
    
    result = runner.run(rho_init)
    assert os.path.exists("convergencia.dat"), "convergencia.dat no fue creado"
    
    # Limpiar
    os.remove("convergencia.dat")
