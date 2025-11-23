# -*- coding: utf-8 -*-
import time
from scf.config import Config
from scf.grid import Grid1D
from scf.potential import ExternalPotentialParams
from scf.scf_loop import InitialDensityParams, InitialDensityBuilder, SCFRunner
from scf.io_utils import IOUtils

J2eV = 1.602176634e-19  # Factor de conversión J → eV

def ask_float(msg):
    return float(input(msg))

def ask_int(msg):
    return int(input(msg))

def main():
    cfg = Config()

    # Entrada interactiva
    k = ask_float("Ingrese el valor de k: ")
    m = ask_float("Ingrese el valor de m: ")
    rho0 = ask_float("Ingrese el valor de rho_0: ")

    print("\nSeleccione el potencial externo:")
    print("1. v(x) = alpha * |x|")
    print("2. v(x) = beta * (x^2 - a^2)^2 (doble pozo)")
    print("3. v(x) = 0.5 * m * omega_ext^2 * x^2 (armónico)")
    print("4. v(x) = -V0 en |x| < L/2 (pozo cuadrado finito)")
    potential_type = ask_int("Opción: ")

    alpha = 1.0
    beta = 1.0
    a_dw = 1.0
    omega_ext = 1.0
    V0 = 1.0
    L = 2.0

    if potential_type == 1:
        alpha = ask_float("Ingrese el valor de alpha (coeficiente de |x| en v(x)): ")
    elif potential_type == 2:
        beta = ask_float("Ingrese el valor de beta (profundidad del doble pozo): ")
        a_dw = ask_float("Ingrese el valor de a (separación entre pozos): ")
    elif potential_type == 3:
        omega_ext = ask_float("Ingrese el valor de omega_ext (frecuencia del oscilador): ")
    elif potential_type == 4:
        V0 = ask_float("Ingrese el valor de V0 (profundidad del pozo cuadrado): ")
        L = ask_float("Ingrese el ancho L del pozo cuadrado: ")

    print("\nSeleccione tipo de densidad inicial:")
    print("1. Gaussiana\n2. Bimodal\n3. Lorentziana")
    guess_type = ask_int("Opción: ")

    # Preparación de malla y densidad inicial
    grid = Grid1D(cfg.N, cfg.x_min, cfg.x_max)
    init = InitialDensityParams(
        guess_type=guess_type,
        Np=cfg.Np,
        sigma=cfg.sigma,
        offset=cfg.offset
    )
    rho_init = InitialDensityBuilder.build(grid, init)

    # Parámetros de potencial externo
    ext_params = ExternalPotentialParams(
        potential_type=potential_type,
        alpha=alpha,
        beta=beta,
        a_dw=a_dw,
        omega_ext=omega_ext,
        V0=V0,
        L=L
    )

    # Ejecutar SCF y medir tiempo
    runner = SCFRunner(cfg, grid, ext_params, m=m, k=k, rho0=rho0)
    start_time = time.time()
    result = runner.run(rho_init)
    end_time = time.time()
    total_time = end_time - start_time

    mu = result["mu"]/J2eV
    E = result["E"]/J2eV
    rho = result["rho"]/J2eV

    # Salidas en consola
    print(f"\nEnergía química (mu): {mu} eV")
    print(f"Energía del estado fundamental (E): {E} eV")
    print(f"Tiempo total de cálculo: {total_time:.2f} segundos.")

    # Archivos y gráficos
    IOUtils.save_density(grid.x, rho, "rho_vs_x.dat")
    IOUtils.write_gnuplot_density_script(E, potential_type, "plot.gp")
    IOUtils.write_gnuplot_convergence_script("convergencia.gp")
    IOUtils.run_gnuplot("plot.gp", "convergencia.gp")
    IOUtils.open_eps("rho_vs_x.eps")
    IOUtils.open_eps("convergencia.eps")

if __name__ == "__main__":
    main()
