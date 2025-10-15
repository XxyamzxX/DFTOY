# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import threading
import queue
import os
from pathlib import Path

# Importa tu motor numérico ya modular
from scf.config import Config
from scf.grid import Grid1D
from scf.potential import ExternalPotentialParams
from scf.scf_loop import InitialDensityParams, InitialDensityBuilder, SCFRunner
from scf.io_utils import IOUtils

from .widgets import ParamsPanel, DensityPanel, PlotPreview, Timeline
from .bridge import run_scf_once

APP_TITLE = "DFT 1D — v1.0"
APP_GEOM = "1100x700"

class DFTApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry(APP_GEOM)
        self.minsize(960, 640)

        # Estado compartido
        self.cfg = Config()
        self.params = {
            "k": 1.0,
            "m": 1.0,
            "rho0": 1.0,
            "potential_type": 3,
            "alpha": 1.0,
            "beta": 1.0,
            "a_dw": 1.0,
            "omega_ext": 1.0,
            "V0": 1.0,
            "L": 2.0,
            "guess_type": 1,
            "sigma": 1.0,
            "offset": 1.0,
        }

        # Cola para mensajes y progreso
        self.msg_queue = queue.Queue()

        # Construcción UI
        self._build_menu()
        self._build_layout()

        # Hilo de trabajo (cuando corresponda)
        self.worker = None
        self.running = False

        # Bindeos
        self.protocol("WM_DELETE_WINDOW", self.on_exit)

    def _build_menu(self):
        menubar = tk.Menu(self)
        # Archivo
        m_file = tk.Menu(menubar, tearoff=0)
        m_file.add_command(label="Nuevo", command=self.on_new)
        m_file.add_command(label="Abrir...", command=self.on_open)
        m_file.add_command(label="Guardar", command=self.on_save)
        m_file.add_separator()
        m_file.add_command(label="Exportar ρ(x)...", command=self.on_export_density)
        m_file.add_separator()
        m_file.add_command(label="Salir", command=self.on_exit)
        menubar.add_cascade(label="Archivo", menu=m_file)
        # Edición
        m_edit = tk.Menu(menubar, tearoff=0)
        m_edit.add_command(label="Preferencias", command=self.on_prefs)
        menubar.add_cascade(label="Edición", menu=m_edit)
        # Simulación
        m_sim = tk.Menu(menubar, tearoff=0)
        m_sim.add_command(label="Ejecutar", command=self.on_run)
        m_sim.add_command(label="Detener", command=self.on_stop)
        menubar.add_cascade(label="Simulación", menu=m_sim)
        # Resultados
        m_res = tk.Menu(menubar, tearoff=0)
        m_res.add_command(label="Potencial externo", command=lambda: self.timeline.log("Mostrar v_ext no implementado"))
        m_res.add_command(label="Potencial efectivo", command=lambda: self.timeline.log("Mostrar v_eff no implementado"))
        m_res.add_command(label="Estados excitados", command=lambda: self.timeline.log("Estados excitados no implementados"))
        m_res.add_command(label="Función de onda", command=lambda: self.open_eps("rho_vs_x.eps"))
        m_res.add_command(label="Convergencia", command=lambda: self.open_eps("convergencia.eps"))
        menubar.add_cascade(label="Resultados", menu=m_res)
        # Ayuda
        m_help = tk.Menu(menubar, tearoff=0)
        m_help.add_command(label="Manual", command=self.on_help)
        menubar.add_cascade(label="Ayuda", menu=m_help)
        self.config(menu=menubar)  # barra de menú estándar Tk [web:64][web:60]

    def _build_layout(self):
        # Contenedor principal con PanedWindow: izquierda (parámetros) / derecha (gráfica)
        main_h = tk.PanedWindow(self, orient=tk.HORIZONTAL, sashrelief=tk.RAISED, showhandle=True)
        main_h.pack(fill=tk.BOTH, expand=True)

        # Panel izquierdo: parámetros y densidad
        left_frame = ttk.Frame(main_h, padding=6)
        left_frame.columnconfigure(0, weight=1)
        left_frame.rowconfigure(0, weight=1)
        left_frame.rowconfigure(1, weight=1)
        self.params_panel = ParamsPanel(left_frame, self.params, on_change=self.on_params_changed)
        self.params_panel.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)
        self.density_panel = DensityPanel(left_frame, self.params, on_change=self.on_params_changed)
        self.density_panel.grid(row=1, column=0, sticky="nsew", padx=4, pady=4)
        main_h.add(left_frame)

        # Panel derecho vertical: gráfico arriba, timeline abajo
        right_v = tk.PanedWindow(main_h, orient=tk.VERTICAL, sashrelief=tk.RAISED, showhandle=True)
        main_h.add(right_v)

        plot_frame = ttk.Frame(right_v, padding=6)
        plot_frame.columnconfigure(0, weight=1)
        plot_frame.rowconfigure(0, weight=1)
        self.plot_preview = PlotPreview(plot_frame)
        self.plot_preview.grid(row=0, column=0, sticky="nsew")
        right_v.add(plot_frame)

        bottom_frame = ttk.Frame(right_v, padding=6)
        bottom_frame.columnconfigure(0, weight=1)
        bottom_frame.rowconfigure(0, weight=1)
        self.timeline = Timeline(bottom_frame)
        self.timeline.grid(row=0, column=0, sticky="nsew")
        right_v.add(bottom_frame)

        # Temporizador para extraer mensajes de la cola y actualizar timeline
        self.after(150, self._poll_msgs)  # patrón común para “status bar / timeline” en Tk [web:76][web:71]

    # Handlers de menú y botones
    def on_new(self):
        self.params_panel.reset_defaults()
        self.density_panel.reset_defaults()
        self.timeline.log("Nuevo proyecto listo.")

    def on_open(self):
        path = filedialog.askopenfilename(title="Abrir parámetros", filetypes=[("Archivos TXT", "*.txt"), ("Todos", "*.*")])
        if not path:
            return
        try:
            with open(path, "r") as f:
                for line in f:
                    if "=" in line:
                        k, v = line.strip().split("=", 1)
                        self.params[k.strip()] = float(v.strip())
            self.params_panel.refresh_from(self.params)
            self.density_panel.refresh_from(self.params)
            self.timeline.log(f"Parámetros cargados desde {path}")
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo abrir: {e}")

    def on_save(self):
        path = filedialog.asksaveasfilename(title="Guardar parámetros", defaultextension=".txt")
        if not path:
            return
        try:
            with open(path, "w") as f:
                for k, v in self.params.items():
                    f.write(f"{k}={v}\n")
            self.timeline.log(f"Parámetros guardados en {path}")
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo guardar: {e}")

    def on_export_density(self):
        target = filedialog.asksaveasfilename(title="Exportar rho_vs_x.dat", defaultextension=".dat")
        if target and Path("rho_vs_x.dat").exists():
            try:
                Path(target).write_text(Path("rho_vs_x.dat").read_text(encoding="utf-8"), encoding="utf-8")
                self.timeline.log(f"ρ(x) exportado a {target}")
            except Exception as e:
                messagebox.showerror("Error", f"No se pudo exportar: {e}")
        else:
            self.timeline.log("No existe rho_vs_x.dat todavía.")

    def on_prefs(self):
        messagebox.showinfo("Preferencias", "No implementado aún.")

    def on_help(self):
        messagebox.showinfo("Manual", "Interfaz DFT 1D.\nUse Simulación > Ejecutar para correr el SCF.")

    def on_params_changed(self, new_params: dict):
        self.params.update(new_params)

    def on_run(self):
        if self.running:
            self.timeline.log("Una simulación ya está en ejecución...")
            return
        self.running = True
        self.timeline.log("Iniciando simulación SCF...")
        # Hilo de trabajo; mantiene UI responsiva
        self.worker = threading.Thread(target=run_scf_once, args=(self.params.copy(), self.msg_queue), daemon=True)
        self.worker.start()

    def on_stop(self):
        self.timeline.log("Detener solicitado: se detendrá al finalizar la iteración actual.")

    def _poll_msgs(self):
        try:
            while True:
                msg = self.msg_queue.get_nowait()
                kind = msg.get("kind", "log")
                if kind == "log":
                    self.timeline.log(msg["text"])
                elif kind == "done":
                    mu = msg.get("mu")
                    E = msg.get("E")
                    t = msg.get("time")
                    self.timeline.log(f"Listo. μ={mu:.6f}, E={E:.6f}, t={t:.3f}s")
                    # Intentar generar/abrir gráficos EPS y previsualizar PNG
                    IOUtils.run_gnuplot("plot.gp", "convergencia.gp")
                    self.plot_preview.refresh_from_eps("rho_vs_x.eps")
                    self.running = False
                elif kind == "error":
                    messagebox.showerror("Error", msg["text"])
                    self.running = False
        except queue.Empty:
            pass
        # rearmar polling
        self.after(150, self._poll_msgs)

    def open_eps(self, path):
        IOUtils.open_eps(path)

    def on_exit(self):
        if self.running:
            if not messagebox.askokcancel("Salir", "Hay una simulación en curso. ¿Salir de todos modos?"):
                return
        self.destroy()

if __name__ == "__main__":
    app = DFTApp()
    app.mainloop()
