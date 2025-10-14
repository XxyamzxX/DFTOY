# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk  # pip install pillow
import subprocess
from pathlib import Path

class LabeledEntry(ttk.Frame):
    def __init__(self, master, label, init_val, width=10, on_change=None):
        super().__init__(master)
        ttk.Label(self, text=label).grid(row=0, column=0, sticky="w")
        self.var = tk.StringVar(value=str(init_val))
        self.entry = ttk.Entry(self, textvariable=self.var, width=width)
        self.entry.grid(row=0, column=1, sticky="e")
        self.on_change = on_change
        self.var.trace_add("write", self._changed)

    def _changed(self, *args):
        if self.on_change:
            self.on_change()

    def value(self, cast=float):
        try:
            return cast(self.var.get())
        except Exception:
            return None

    def set_value(self, v):
        self.var.set(str(v))

class ParamsPanel(ttk.Labelframe):
    def __init__(self, master, params, on_change=None):
        super().__init__(master, text="Agregar parámetros iniciales")
        self.on_change = on_change
        self._build(params)

    def _build(self, params):
        # Filas de parámetros físicos
        self.inputs = {}
        grid = [
            ("k", params["k"]),
            ("m", params["m"]),
            ("rho0", params["rho0"]),
            ("omega_ext", params["omega_ext"]),
            ("alpha", params["alpha"]),
            ("beta", params["beta"]),
            ("a_dw", params["a_dw"]),
            ("V0", params["V0"]),
            ("L", params["L"]),
            ("potential_type", params["potential_type"]),
        ]
        for r, (name, val) in enumerate(grid):
            w = LabeledEntry(self, name, val, on_change=self._notify)
            w.grid(row=r, column=0, sticky="ew", padx=4, pady=2)
            self.inputs[name] = w
        self.columnconfigure(0, weight=1)

    def _notify(self):
        if self.on_change:
            self.on_change(self.values())

    def values(self):
        vals = {k: self.inputs[k].value(float) for k in self.inputs}
        vals["potential_type"] = int(self.inputs["potential_type"].value(int))
        return vals

    def reset_defaults(self):
        for k, w in self.inputs.items():
            w.set_value(w.value())

    def refresh_from(self, params):
        for k, w in self.inputs.items():
            w.set_value(params[k])

class DensityPanel(ttk.Labelframe):
    def __init__(self, master, params, on_change=None):
        super().__init__(master, text="Agregar función de densidad")
        self.on_change = on_change
        self._build(params)

    def _build(self, params):
        self.inputs = {}
        grid = [
            ("guess_type", params["guess_type"]),
            ("sigma", params["sigma"]),
            ("offset", params["offset"]),
        ]
        for r, (name, val) in enumerate(grid):
            w = LabeledEntry(self, name, val, on_change=self._notify)
            w.grid(row=r, column=0, sticky="ew", padx=4, pady=2)
            self.inputs[name] = w
        self.columnconfigure(0, weight=1)

    def _notify(self):
        if self.on_change:
            data = {k: self.inputs[k].value(float) for k in self.inputs}
            data["guess_type"] = int(self.inputs["guess_type"].value(int))
            self.on_change(data)

    def reset_defaults(self):
        for k, w in self.inputs.items():
            w.set_value(w.value())

    def refresh_from(self, params):
        for k, w in self.inputs.items():
            if k in params:
                w.set_value(params[k])

class PlotPreview(ttk.Labelframe):
    def __init__(self, master):
        super().__init__(master, text="Gráfica del DFT")
        self.canvas = tk.Label(self)
        self.canvas.pack(fill=tk.BOTH, expand=True)
        self._img = None

    def refresh_from_eps(self, eps_path):
        # Convertir EPS -> PNG temporal (requiere ghostscript o imagemagick)
        eps = Path(eps_path)
        if not eps.exists():
            return
        png = eps.with_suffix(".png")
        # Intentar usar ghostscript vía ImageMagick 'convert' si está disponible
        try:
            subprocess.run(["convert", str(eps), str(png)], check=True)
        except Exception:
            # Si no hay convert, no hacemos preview pero no rompemos
            return
        if png.exists():
            img = Image.open(png)
            img.thumbnail((800, 500))
            self._img = ImageTk.PhotoImage(img)
            self.canvas.configure(image=self._img)

class Timeline(ttk.Labelframe):
    def __init__(self, master):
        super().__init__(master, text="Línea de tiempo de mensajes y alertas")
        self.text = tk.Text(self, height=7, wrap="word")
        self.text.pack(fill=tk.BOTH, expand=True)
        self.text.insert("end", "Estado de ejecución...\n")
        self.text.config(state="disabled")

    def log(self, msg: str):
        self.text.config(state="normal")
        self.text.insert("end", msg + "\n")
        self.text.see("end")
        self.text.config(state="disabled")
