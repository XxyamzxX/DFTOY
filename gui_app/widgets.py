# gui_app/widgets.py
# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class LabeledEntry(ttk.Frame):
    def __init__(self, master, label, init_val, width=10, on_change=None, anchor="w"):
        super().__init__(master)
        self.label = ttk.Label(self, text=label, anchor=anchor, width=24)
        self.label.grid(row=0, column=0, sticky="ew", padx=(0,6))
        self.var = tk.StringVar(value=str(init_val))
        self.entry = ttk.Entry(self, textvariable=self.var, width=12)
        self.entry.grid(row=0, column=1, sticky="ew")
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.on_change = on_change
        if on_change:
            self.var.trace_add("write", lambda *_: on_change())

    def value(self, cast=float):
        try:
            return cast(self.var.get())
        except Exception:
            return None

    def set_value(self, v):
        self.var.set(str(v))

class PotentialCombo(ttk.Frame):
    NAMES = ["Lineal |x|", "Doble pozo", "Armónico", "Pozo cuadrado"]
    IDS   = [1, 2, 3, 4]
    def __init__(self, master, init_idx=2, on_change=None):
        super().__init__(master)
        ttk.Label(self, text="Potencial").grid(row=0, column=0, sticky="w", padx=(0,6))
        self.var = tk.StringVar(value=self.NAMES[init_idx])
        self.combo = ttk.Combobox(self, textvariable=self.var, values=self.NAMES, state="readonly", width=20)
        self.combo.grid(row=0, column=1, sticky="ew")
        if on_change:
            self.combo.bind("<<ComboboxSelected>>", lambda e: on_change())
        self.columnconfigure(1, weight=1)

    def get_id(self):
        return self.IDS[self.NAMES.index(self.var.get())]

    def set_by_id(self, pid: int):
        if pid in self.IDS:
            self.var.set(self.NAMES[self.IDS.index(pid)])

class DensityCombo(ttk.Frame):
    NAMES = ["Gaussiana", "Bimodal", "Lorentziana"]
    IDS   = [1, 2, 3]
    def __init__(self, master, init_idx=0, on_change=None):
        super().__init__(master)
        ttk.Label(self, text="Tipo de densidad").grid(row=0, column=0, sticky="w", padx=(0,6))
        self.var = tk.StringVar(value=self.NAMES[init_idx])
        self.combo = ttk.Combobox(self, textvariable=self.var, values=self.NAMES, state="readonly", width=20)
        self.combo.grid(row=0, column=1, sticky="ew")
        if on_change:
            self.combo.bind("<<ComboboxSelected>>", lambda e: on_change())
        self.columnconfigure(1, weight=1)
    def get_id(self):
        return self.IDS[self.NAMES.index(self.var.get())]
    def set_by_id(self, gid: int):
        if gid in self.IDS:
            self.var.set(self.NAMES[self.IDS.index(gid)])

class PlotArea(ttk.Labelframe):
    def __init__(self, master, title="Gráfica del DFT"):
        super().__init__(master, text=title)
        self.fig = Figure(figsize=(6,3.6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.grid(True, linestyle="--", alpha=0.3)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def clear(self):
        self.ax.cla()
        self.ax.grid(True, linestyle="--", alpha=0.3)

    def plot_xy(self, x, y, label=None, xlabel="x", ylabel="", title=""):
        self.clear()
        self.ax.plot(x, y, lw=1.8, label=label)
        if label:
            self.ax.legend(loc="best")
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.set_title(title)
        self.canvas.draw_idle()

class Timeline(ttk.Labelframe):
    def __init__(self, master, title="Línea de tiempo de mensajes y alertas"):
        super().__init__(master, text=title)
        self.text = tk.Text(self, height=8, wrap="word")
        self.text.pack(fill=tk.BOTH, expand=True)
        self.text.insert("end", "Estado de ejecución...\n")
        self.text.config(state="disabled")

    def log(self, msg: str):
        self.text.config(state="normal")
        self.text.insert("end", msg + "\n")
        self.text.see("end")
        self.text.config(state="disabled")

    def clear(self):
        self.text.config(state="normal")
        self.text.delete("1.0", "end")
        self.text.config(state="disabled")
