"""
UI.py
=====
Small Tkinter front-end for a SINGLE insertion.

Pick a .tif volume, choose one of the probes from config.py (or type your own
dimensions), and run one insertion at a chosen (x, y). It reports the
intersected vessel-voxel count and opens the 3-D view.

This is an exploration tool; the published numbers come from the notebooks.

    python UI.py

All dimensions are diameters in micrometres, as everywhere else.
"""

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from config import ELECTRODES, DEPTH_LIMIT, FIG_CROP_X, FIG_CROP_Y
from model_3D_visualization import visualize_cone_pyvista, cropping_img
from Volume_bleeding import load_volume, simulate_cone_insertion, find_first_black_pixel_slice

# Presets come straight from config.py, so the GUI can never disagree with the pipeline.
PRESETS = {"Custom (Manual)": {}, **{cfg["label"]: cfg for cfg in ELECTRODES.values()}}

# Field label -> config key, in the order shown in the window.
FIELDS = [
    ("Shank Length",        "shank_length"),
    ("Shank Base Diameter", "shank_base_d"),
    ("Shank Top Diameter",  "shank_top_d"),
    ("Tip Length",          "tip_length"),
    ("Tip Base Diameter",   "tip_base_d"),
    ("Tip Top Diameter",    "tip_top_d"),
]


def run_gui():
    """Build and start the Tkinter window. Blocks until the window is closed."""

    def browse_file():
        path = filedialog.askopenfilename(filetypes=[("TIFF files", "*.tif *.tiff")])
        entry_file.delete(0, tk.END)
        entry_file.insert(0, path)

    def apply_preset(event=None):
        config = PRESETS.get(probe_combo.get(), {})
        for (_, key), entry in zip(FIELDS, entries):
            if key in config:
                entry.delete(0, tk.END)
                entry.insert(0, str(config[key]))

    def run_simulation():
        try:
            img = load_volume(entry_file.get())
            g = {key: float(entry.get()) for (_, key), entry in zip(FIELDS, entries)}
            x, y = int(entry_x.get()), int(entry_y.get())
            geom = dict(shank_length=g["shank_length"],
                        shank_base_diameter=g["shank_base_d"], shank_top_diameter=g["shank_top_d"],
                        tip_length=g["tip_length"],
                        tip_base_diameter=g["tip_base_d"], tip_top_diameter=g["tip_top_d"])

            start = find_first_black_pixel_slice(img, x, y)
            result = simulate_cone_insertion(img, x, y, **geom, start_slice=start, depth_limit=DEPTH_LIMIT)

            slab, cx, cy = cropping_img(img, x, y, FIG_CROP_X, FIG_CROP_Y, start)
            visualize_cone_pyvista(slab, cx, cy, **geom, start_slice=start, depth_limit=DEPTH_LIMIT, ui=1)

            messagebox.showinfo("Result", f"Total vessel voxels intersected: {int(result)}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    window = tk.Tk()
    window.title("Electrode Insertion Simulation")

    tk.Label(window, text="Load File").grid(row=0, column=0)
    entry_file = tk.Entry(window, width=40)
    entry_file.grid(row=0, column=1)
    tk.Button(window, text="Browse", command=browse_file).grid(row=0, column=2)

    tk.Label(window, text="Preset Electrode").grid(row=1, column=0)
    probe_combo = ttk.Combobox(window, values=list(PRESETS), state="readonly")
    probe_combo.current(0)
    probe_combo.grid(row=1, column=1)
    probe_combo.bind("<<ComboboxSelected>>", apply_preset)

    entries = []
    for i, (label, _) in enumerate(FIELDS):
        tk.Label(window, text=f"{label} (um)").grid(row=i + 2, column=0)
        entry = tk.Entry(window)
        entry.grid(row=i + 2, column=1)
        entries.append(entry)

    row = len(FIELDS) + 2
    tk.Label(window, text="X Center").grid(row=row, column=0)
    entry_x = tk.Entry(window)
    entry_x.grid(row=row, column=1)
    tk.Label(window, text="Y Center").grid(row=row + 1, column=0)
    entry_y = tk.Entry(window)
    entry_y.grid(row=row + 1, column=1)

    tk.Button(window, text="Run Simulation", command=run_simulation).grid(row=row + 3, column=1)
    window.mainloop()


if __name__ == "__main__":
    run_gui()
