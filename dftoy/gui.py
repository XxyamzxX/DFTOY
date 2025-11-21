# gui_app/gui.py
# -*- coding: utf-8 -*-

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import importlib.resources as resources
import threading, queue, json, subprocess, os, sys
import numpy as np
from pathlib import Path

from scf.config import Config
from scf.grid import Grid1D
from scf.potential import ExternalPotentialParams, Potentials
from scf.solver import Solver
from .bridge import run_scf_once
from .widgets import LabeledEntry, PotentialCombo, DensityCombo, PlotArea, Timeline

APP_TITLE = "DFT 1D — v1.3"
APP_GEOM  = "1100x700"
EV_PER_J  = 1.0 / 1.602176634e-19  # J -> eV


class DFTApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry(APP_GEOM)
        self.minsize(960, 640)

        # Estado
        self.cfg = Config()
        self.params = {
            "k": 1.0, "m": 1.0, "rho0": 1.0,
            "potential_type": 3, "alpha": 1.0, "beta": 1.0, "a_dw": 1.0,
            "omega_ext": 1.0, "V0": 1.0, "L": 2.0,
            "guess_type": 1, "sigma": 1.0, "offset": 1.0,
        }
        self.current_file: Path | None = None
        self.msg_queue: queue.Queue = queue.Queue()
        self.worker: threading.Thread | None = None
        self.cancel_evt: threading.Event | None = None
        self.running: bool = False
        self.has_results: bool = False
        self.current_view: str | None = None  # "vext","veff","states","density","conv"

        # UI
        self._build_menu()
        self._build_layout()
        self.protocol("WM_DELETE_WINDOW", self.on_exit)
        self.after(150, self._poll_queue)

    # ================= Unidades (eV) =================
    def _energy_scale_to_ev(self) -> float:
        """
        Heurística: si hbar ~ 1.054e-34 (J·s) convierte J->eV.
        Si hbar ~ 6.582e-16 (eV·s) deja igual. Por defecto: 1.
        """
        hb = float(self.cfg.hbar)
        if 1e-35 < hb < 1e-33:   # SI (J·s)
            return EV_PER_J
        if 1e-15 < hb < 1e-14:   # eV·s
            return 1.0
        return 1.0

    def _to_ev(self, y):
        return np.asarray(y) * self._energy_scale_to_ev()

    # ================= Menús =================
    def _build_menu(self):
        menubar = tk.Menu(self)

        # Archivo
        m_file = tk.Menu(menubar, tearoff=0)
        m_file.add_command(label="Nuevo",   command=self.on_new)
        m_file.add_command(label="Abrir...", command=self.on_open)
        m_file.add_command(label="Guardar",  command=self.on_save)
        m_file.add_separator()
        m_file.add_command(label="Exportar ρ(x)...", command=self.on_export_density)
        m_file.add_separator()
        m_file.add_command(label="Salir", command=self.on_exit)
        menubar.add_cascade(label="Archivo", menu=m_file)

        # Edición
        m_edit = tk.Menu(menubar, tearoff=0)
        m_edit.add_command(label="Limpiar terminal", command=self.clear_console)
        menubar.add_cascade(label="Edición", menu=m_edit)

        # Simulación
        m_sim = tk.Menu(menubar, tearoff=0)
        m_sim.add_command(label="Ejecutar SCF", command=self.start_scf)
        m_sim.add_command(label="Cancelar",     command=self.stop_scf)
        menubar.add_cascade(label="Simulación", menu=m_sim)

        # Resultados
        self.m_res = tk.Menu(menubar, tearoff=0)
        self.idx_vext  = 0; self.m_res.add_command(label="Potencial externo",  command=self.show_v_ext)
        self.idx_veff  = 1; self.m_res.add_command(label="Potencial efectivo", command=self.show_v_eff)
        self.idx_state = 2; self.m_res.add_command(label="Estados (n=0..3)",   command=lambda: self.show_states(4))
        self.idx_rho   = 3; self.m_res.add_command(label="Densidad ρ(x)",      command=self.show_density)
        self.idx_conv  = 4; self.m_res.add_command(label="Convergencia",       command=self.show_convergence)
        menubar.add_cascade(label="Resultados", menu=self.m_res)

        # Ayuda
        m_help = tk.Menu(menubar, tearoff=0)
        m_help.add_command(label="Manual de usuario", command=self.open_manual)
        m_help.add_separator()
        m_help.add_command(label="Acerca de…", command=lambda: messagebox.showinfo("Acerca de", APP_TITLE))
        menubar.add_cascade(label="Ayuda", menu=m_help)

        self.config(menu=menubar)
        self._set_results_enabled(False)

    def _set_results_enabled(self, enabled: bool):
        state = "normal" if enabled else "disabled"
        for i in [self.idx_vext, self.idx_veff, self.idx_state, self.idx_rho, self.idx_conv]:
            self.m_res.entryconfig(i, state=state)

    # ================= Layout =================
    def _build_layout(self):
        main_h = tk.PanedWindow(self, orient=tk.HORIZONTAL, sashrelief=tk.RAISED, showhandle=True)
        main_h.pack(fill=tk.BOTH, expand=True)

        # Panel izquierdo
        left = ttk.Labelframe(main_h, text="Parámetros")
        self.fields: dict[str, LabeledEntry] = {}

        def add_row(key, label):
            w = LabeledEntry(left, label, self.params[key], on_change=None)
            w.grid(sticky="ew", padx=6, pady=3)
            self.fields[key] = w

        add_row("k",         "k (constante elástica)")
        add_row("m",         "m (masa)")
        add_row("rho0",      "ρ₀ (densidad ref.)")
        add_row("omega_ext", "ω_ext (frecuencia)")
        add_row("alpha",     "α (lineal)")
        add_row("beta",      "β (doble pozo)")
        add_row("a_dw",      "a_dw (ancho del pozo)")
        add_row("V0",        "V0 (profundidad)")
        add_row("L",         "L (ancho caja)")

        ttk.Label(left, text="Tipo de potencial").grid(sticky="w", padx=6, pady=(8,0))
        self.pcombo = PotentialCombo(left, init_idx=self.params["potential_type"]-1, on_change=None)
        self.pcombo.grid(sticky="ew", padx=6, pady=3)

        dens = ttk.Labelframe(left, text="Densidad inicial")
        dens.grid(sticky="ew", padx=6, pady=(8,6))
        self.d_combo  = DensityCombo(dens, init_idx=self.params["guess_type"]-1)
        self.d_sigma  = LabeledEntry(dens, "σ (ancho)",  self.params["sigma"])
        self.d_offset = LabeledEntry(dens, "offset",     self.params["offset"])
        self.d_combo.grid(sticky="ew", padx=4, pady=2)
        self.d_sigma.grid(sticky="ew", padx=4, pady=2)
        self.d_offset.grid(sticky="ew", padx=4, pady=2)

        left.columnconfigure(0, weight=1)
        main_h.add(left, stretch="always")

        # Panel derecho
        right_v = tk.PanedWindow(main_h, orient=tk.VERTICAL, sashrelief=tk.RAISED, showhandle=True)
        self.plot = PlotArea(right_v, title="Gráfica")
        right_v.add(self.plot, stretch="always")
        self.timeline = Timeline(right_v)
        right_v.add(self.timeline, stretch="always")
        main_h.add(right_v, stretch="always")

    # ================= Helpers =================
    def _current_grid(self) -> Grid1D:
        return Grid1D(self.cfg.N, self.cfg.x_min, self.cfg.x_max)

    def _read_fields_to_params(self):
        def get_num(widget, cast=float, fallback=None):
            try:
                return cast(widget.value(cast))
            except Exception:
                return fallback

        for key in ["k","m","rho0","omega_ext","alpha","beta","a_dw","V0","L"]:
            self.params[key] = get_num(self.fields[key], float, self.params[key])
        self.params["potential_type"] = int(self.pcombo.get_id())
        self.params["guess_type"]     = int(self.d_combo.get_id())
        self.params["sigma"]  = float(get_num(self.d_sigma,  float, self.params["sigma"]))
        self.params["offset"] = float(get_num(self.d_offset, float, self.params["offset"]))

    def _refresh_fields_from_params(self):
        for k, w in self.fields.items():
            w.set_value(self.params[k])
        self.pcombo.set_by_id(int(self.params["potential_type"]))
        self.d_combo.set_by_id(int(self.params["guess_type"]))
        self.d_sigma.set_value(self.params["sigma"])
        self.d_offset.set_value(self.params["offset"])

    # ================= Archivo =================
    def ensure_dft_ext(self, path_str: str) -> str:
        p = Path(path_str)
        if p.suffix.lower() != ".dft":
            p = p.with_suffix(".dft")
        return str(p)

    def on_new(self):
        self.params = {
            "k": 1.0, "m": 1.0, "rho0": 1.0,
            "potential_type": 3, "alpha": 1.0, "beta": 1.0, "a_dw": 1.0,
            "omega_ext": 1.0, "V0": 1.0, "L": 2.0,
            "guess_type": 1, "sigma": 1.0, "offset": 1.0,
        }
        self.current_file = None
        self._refresh_fields_from_params()
        self.plot.clear(); self.plot.canvas.draw_idle()
        self.timeline.clear(); self.timeline.log("Nuevo proyecto listo.")
        self.has_results = False
        self._set_results_enabled(False)

    def on_open(self):
        path = filedialog.askopenfilename(
            title="Abrir proyecto",
            filetypes=[("DFT (*.dft)", "*.dft"), ("Todos", "*.*")]
        )
        if not path:
            return
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
            if isinstance(data, dict):
                self.params.update(data)
                self._refresh_fields_from_params()
                self.current_file = Path(path)
                self.timeline.log(f"Proyecto cargado: {path}")
                self.has_results = False
                self._set_results_enabled(False)
            else:
                messagebox.showerror("Error", "Archivo .dft inválido.")
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo abrir: {e}")

    def on_save(self):
        self._read_fields_to_params()
        if not self.current_file:
            path = filedialog.asksaveasfilename(
                title="Guardar proyecto",
                defaultextension=".dft",
                filetypes=[("DFT (*.dft)", "*.dft"), ("Todos", "*.*")],
                initialfile="proyecto.dft",
            )
            if not path:
                return
            path = self.ensure_dft_ext(path)
            self.current_file = Path(path)
        try:
            self.current_file.write_text(json.dumps(self.params, indent=2), encoding="utf-8")
            self.timeline.log(f"Proyecto guardado en {self.current_file}")
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo guardar: {e}")

    def on_export_density(self):
        src = Path("rho_vs_x.dat")
        if not src.exists():
            self.timeline.log("Aún no existe rho_vs_x.dat; ejecute la simulación.")
            return
        target = filedialog.asksaveasfilename(
            title="Exportar ρ(x)",
            defaultextension=".dat",
            filetypes=[("Datos (*.dat)", "*.dat"), ("CSV (*.csv)", "*.csv"), ("Todos", "*.*")],
            initialfile="rho_vs_x.dat",
        )
        if not target:
            return
        try:
            Path(target).write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
            self.timeline.log(f"ρ(x) exportado a {target}")
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo exportar: {e}")

    # ================= Simulación =================
    def start_scf(self):
        if self.running:
            self.timeline.log("Ya hay una simulación en ejecución.")
            return
        self._read_fields_to_params()
        self.running = True
        self.has_results = False
        self._set_results_enabled(False)
        self.cancel_evt = threading.Event()
        self.timeline.log("Iniciando SCF, preparando datos de malla...")
        self.worker = threading.Thread(
            target=run_scf_once,
            args=(self.params.copy(), self.msg_queue, self.cancel_evt),
            daemon=True
        )
        self.worker.start()

    def stop_scf(self):
        if not self.running:
            self.timeline.log("No hay simulación en curso.")
            return
        if self.cancel_evt:
            self.cancel_evt.set()
        if self.worker:
            self.worker.join(timeout=0.5)
        self.running = False
        self.has_results = False
        self._set_results_enabled(False)
        self.timeline.log("Cancelado.")

    def _poll_queue(self):
        try:
            while True:
                msg = self.msg_queue.get_nowait()
                kind = msg.get("type", "log")
                if kind == "log":
                    self.timeline.log(msg.get("text", ""))
                elif kind == "canceled":
                    self.timeline.log("Simulación cancelada por el usuario.")
                    self.running = False
                    self.has_results = False
                    self._set_results_enabled(False)
                elif kind == "done":
                    mu = msg.get("mu"); E = msg.get("E"); t = msg.get("time")
                    mu_ev = float(self._to_ev(mu)) if mu is not None else None
                    E_ev  = float(self._to_ev(E))  if E  is not None else None
                    self.timeline.log("Simulación finalizada; ya puede ver los resultados en ‘Resultados’.")
                    if mu_ev is not None and E_ev is not None and t is not None:
                        self.timeline.log(f"μ={mu_ev:.6f} eV, E={E_ev:.6f} eV, t={t:.3f} s")
                    self.running = False
                    self.has_results = True
                    self._set_results_enabled(True)
                    self._refresh_current_view()
                elif kind == "error":
                    self.timeline.log(f"Error: {msg.get('text','')}")
                    self.running = False
                    self.has_results = False
                    self._set_results_enabled(False)
        except queue.Empty:
            pass
        self.after(150, self._poll_queue)

    def _refresh_current_view(self):
        if self.current_view == "vext":
            self.show_v_ext()
        elif self.current_view == "veff":
            self.show_v_eff()
        elif self.current_view == "states":
            self.show_states(4)
        elif self.current_view == "conv":
            self.show_convergence()
        else:
            # Por defecto, mostrar densidad al finalizar
            self.show_density()

    def on_exit(self):
        if self.running:
            if not messagebox.askyesno("Salir", "Hay una simulación en curso. ¿Desea salir igualmente?"):
                return
            self.stop_scf()
        self.destroy()

    def clear_console(self):
        self.timeline.clear()
        try:
            if sys.stdout.isatty():
                subprocess.call("cls" if os.name == "nt" else "clear", shell=True)
        except Exception:
            pass
        self.timeline.log("Terminal limpia.")

    # ================= Resultados =================
    def _params_objects(self):
        p = ExternalPotentialParams(
            potential_type=int(self.params["potential_type"]),
            alpha=float(self.params["alpha"]), beta=float(self.params["beta"]),
            a_dw=float(self.params["a_dw"]), omega_ext=float(self.params["omega_ext"]),
            V0=float(self.params["V0"]), L=float(self.params["L"]),
        )
        m = float(self.params["m"]); k = float(self.params["k"]); rho0 = float(self.params["rho0"])
        return p, m, k, rho0

    def show_v_ext(self):
        grid = self._current_grid(); self._read_fields_to_params()
        p, m, _, _ = self._params_objects()
        v_ev = self._to_ev(Potentials.v_ext(grid.x, p, m))
        self.plot.plot_xy(grid.x, v_ev, xlabel="x", ylabel="V_ext(x) [eV]", title="Potencial externo")
        self.current_view = "vext"

    def show_v_eff(self):
        grid = self._current_grid(); self._read_fields_to_params()
        p, m, k, rho0 = self._params_objects()
        omega0, omega = Potentials.compute_omegas(self.cfg.hbar, k, m, rho0)
        vext = Potentials.v_ext(grid.x, p, m)
        veff = Potentials.v_eff(grid.x, vext, self.cfg.hbar, m, omega0, omega)
        self.plot.plot_xy(grid.x, self._to_ev(veff), xlabel="x", ylabel="V_eff(x) [eV]", title="Potencial efectivo")
        self.current_view = "veff"

    def show_states(self, n=4):
        grid = self._current_grid(); self._read_fields_to_params()
        p, m, k, rho0 = self._params_objects()
        omega0, omega = Potentials.compute_omegas(self.cfg.hbar, k, m, rho0)
        vext = Potentials.v_ext(grid.x, p, m)
        veff = Potentials.v_eff(grid.x, vext, self.cfg.hbar, m, omega0, omega)
        dx = grid.dx
        t = Config.kinetic_coeff(self.cfg.hbar, m, dx)
        diag = veff + (-2.0 * t)
        off  = (t * np.ones(grid.N - 1))
        w, V = Solver.eigenpairs_tridiagonal(diag, off, n)

        self.plot.clear()
        scale = 1.0 / (abs(V).max() + 1e-12)
        wev = self._to_ev(w); veff_ev = self._to_ev(veff)
        for i in range(min(n, V.shape[1])):
            self.plot.ax.plot(grid.x, V[:, i] * scale + wev[i], lw=1.2, label=f"n={i}")
        self.plot.ax.plot(grid.x, veff_ev, 'k--', lw=1.0, label="V_eff")
        self.plot.ax.set_xlabel("x"); self.plot.ax.set_ylabel("E [eV]")
        self.plot.ax.set_title("Estados propios (desplazados)")
        self.plot.ax.legend(loc="best"); self.plot.canvas.draw_idle()
        self.current_view = "states"

    def show_density(self):
        if not self.has_results:
            self.timeline.log("No hay resultados de esta sesión; ejecute SCF primero.")
            return
        try:
            xs, ys = [], []
            with open("rho_vs_x.dat","r") as f:
                for line in f:
                    x,y = map(float, line.split()); xs.append(x); ys.append(y)
            self.plot.plot_xy(xs, ys, xlabel="x", ylabel="ρ(x)", title="Densidad electrónica")
            self.current_view = "density"
        except Exception:
            self.timeline.log("No se pudo leer rho_vs_x.dat.")

    def show_convergence(self):
        if not self.has_results:
            self.timeline.log("No hay resultados de esta sesión; ejecute SCF primero.")
            return
        try:
            its, errs = [], []
            with open("convergencia.dat","r") as f:
                for line in f:
                    it, err = map(float, line.split()); its.append(it); errs.append(err)
            self.plot.plot_xy(its, errs, xlabel="Iteración", ylabel="‖Δρ‖", title="Convergencia")
            self.current_view = "conv"
        except Exception:
            self.timeline.log("No se pudo leer convergencia.dat.")

    # ================= Ayuda =================
    def open_manual(self):
        """
        Abre el PDF del manual de usuario empaquetado en:
        dftoy/Manual_Usuario/Manual_usuario.pdf
        Funciona tanto en entorno local como paquete instalado.
        """
        try:
            # 'dftoy.Manual_Usuario' es la carpeta dentro del paquete
            pdf_file = resources.files("dftoy") / "Manual_Usuario" / "Manual_usuario.pdf"

            # Convierte a Path para poder abrirlo
            pdf_path = Path(pdf_file)
            if not pdf_path.exists():
                raise FileNotFoundError(f"No se encontró {pdf_path}")

            # Abrir con el visor predeterminado
            self._open_path(pdf_path)
            if hasattr(self, "timeline"):
                self.timeline.log(f"Abriendo manual de usuario: {pdf_path}")

        except FileNotFoundError:
            messagebox.showerror(
                "Manual no encontrado",
                "No se encontró Manual_usuario.pdf. Verifique que el paquete se instaló correctamente."
            )
        except Exception as e:
            messagebox.showerror(
                "Error al abrir manual",
                f"No se pudo abrir el manual: {e}"
            )



    def _open_path(self, path: Path):
        """
        Abre un archivo con el visor predeterminado del sistema.
        - Windows: os.startfile
        - macOS: open
        - Linux: xdg-open
        - WSL (si está disponible): wslview
        """
        try:
            # WSL: preferir wslview si existe (mejor manejo de rutas /mnt/c)
            if Path("/usr/bin/wslview").exists():
                subprocess.Popen(["wslview", str(path)])
                return
        except Exception:
            pass

        try:
            if sys.platform.startswith("win"):
                os.startfile(str(path))  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                subprocess.Popen(["open", str(path)])
            else:
                subprocess.Popen(["xdg-open", str(path)])
        except Exception as e:
            messagebox.showerror("Error al abrir manual", f"No se pudo abrir el manual:\n{e}")


def main():
    app = DFTApp()
    # Soporte para abrir un .dft con doble clic (python -m gui_app.gui proyecto.dft)
    if len(sys.argv) >= 2:
        try:
            data = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
            if isinstance(data, dict):
                app.params.update(data)
                app._refresh_fields_from_params()
                app.current_file = Path(sys.argv[1])
                app.timeline.log(f"Proyecto cargado: {sys.argv[1]}")
        except Exception as e:
            app.timeline.log(f"No se pudo abrir {sys.argv[1]}: {e}")
    app.mainloop()

if __name__ == "__main__":
    main()
