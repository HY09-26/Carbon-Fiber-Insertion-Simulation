"""
UI.py
=====
Small Tkinter front-end for a SINGLE insertion.

Pick a .tif volume, choose one of the modelled probes (or type your own
dimensions), and run one insertion at a chosen (x, y). It reports the
intersected vessel-voxel count and opens the 3-D view.

This is an exploration tool. The numbers that go into the paper come from
`All_Simulation.ipynb`, which sweeps 100 insertion sites across every animal.

Run it with:

    python UI.py

All dimensions are DIAMETERS in micrometres, matching the specification table
in README.md and the simulation functions, so nothing is converted anywhere.
"""

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import tifffile as tiff
import numpy as np

from model_3D_visualization import visualize_cone_pyvista, cropping_img
from Volume_bleeding import simulate_cone_insertion, find_first_black_pixel_slice

# Probe presets, in micrometres. "shank" is the segment nearer the cortical
# surface, "tip" the deepest one; each tapers linearly from base to top.
# Shank + tip length is 1000 um for every probe.
PRESETS = {
    "Custom (Manual)": {},
    "Carbon Fiber": dict(shank_length=840,    shank_base_d=8.4, shank_top_d=8.4,
                         tip_length=160,      tip_base_d=6.8,   tip_top_d=0.0),
    "Microprobes FMA": dict(shank_length=971.91, shank_base_d=25.0, shank_top_d=25.0,
                            tip_length=28.09,    tip_base_d=25.0,   tip_top_d=0.0),
    "Shuttle": dict(shank_length=990,        shank_base_d=25.0, shank_top_d=25.0,
                    tip_length=10,           tip_base_d=25.0,   tip_top_d=25.0),
    "Blackrock UEA": dict(shank_length=950,  shank_base_d=90.0, shank_top_d=28.0,
                          tip_length=50,     tip_base_d=28.0,   tip_top_d=3.0),
}

# Field label -> preset key. Order here is the order shown in the window.
FIELDS = [
    ("Shank Length",        "shank_length"),
    ("Shank Base Diameter", "shank_base_d"),
    ("Shank Top Diameter",  "shank_top_d"),
    ("Tip Length",          "tip_length"),
    ("Tip Base Diameter",   "tip_base_d"),
    ("Tip Top Diameter",    "tip_top_d"),
]

CROP_MARGIN_X = 80        # half-width of the rendered slab along x, in voxels
CROP_MARGIN_Y = 50        # half-width along y; 50 fills the 100-voxel slab
DEPTH_LIMIT = 1000        # insertion depth, in voxels


def run_gui():
    """Build and start the Tkinter window. Blocks until the window is closed."""
    print("UI launching...")

    def browse_file():
        """File picker for the .tif volume."""
        path = filedialog.askopenfilename(filetypes=[("TIFF files", "*.tif *.tiff")])
        entry_file.delete(0, tk.END)
        entry_file.insert(0, path)

    def apply_preset(event=None):
        """Fill the dimension fields from the selected probe preset."""
        config = PRESETS.get(probe_combo.get(), {})
        for (_, key), entry in zip(FIELDS, entries):
            if key in config:
                entry.delete(0, tk.END)
                entry.insert(0, str(config[key]))

    def run_simulation():
        """Load the volume, run one insertion, report the count, show the 3-D view."""
        try:
            # (z, y, x) as stored on disk -> (z, x, y) as the simulation expects.
            img_data = tiff.imread(entry_file.get())
            img_data = np.transpose(img_data, axes=(0, 2, 1)).astype(np.uint16)
            print(img_data.shape)

            geom = {key: float(entry.get()) for (_, key), entry in zip(FIELDS, entries)}
            x_center = int(entry_x.get())
            y_center = int(entry_y.get())

            start_slice = find_first_black_pixel_slice(img_data, x_center, y_center)

            result = simulate_cone_insertion(
                img_data, x_center, y_center,
                geom["shank_length"], geom["shank_base_d"], geom["shank_top_d"],
                geom["tip_length"], geom["tip_base_d"], geom["tip_top_d"],
                start_slice, DEPTH_LIMIT)

            # Crop to a slab around the insertion so the renderer stays responsive.
            cropped, cx, cy = cropping_img(img_data, x_center, y_center,
                                           CROP_MARGIN_X, CROP_MARGIN_Y, start_slice)

            visualize_cone_pyvista(
                cropped, cx, cy,
                shank_length=geom["shank_length"],
                shank_base_diameter=geom["shank_base_d"],
                shank_top_diameter=geom["shank_top_d"],
                tip_length=geom["tip_length"],
                tip_base_diameter=geom["tip_base_d"],
                tip_top_diameter=geom["tip_top_d"],
                start_slice=start_slice, depth_limit=DEPTH_LIMIT, ui=1)

            messagebox.showinfo("Result", f"Total vessel voxels intersected: {int(result)}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    window = tk.Tk()
    window.title("Electrode Insertion Simulation")

    # File input
    tk.Label(window, text="Load File").grid(row=0, column=0)
    entry_file = tk.Entry(window, width=40)
    entry_file.grid(row=0, column=1)
    tk.Button(window, text="Browse", command=browse_file).grid(row=0, column=2)

    # Probe preset dropdown
    tk.Label(window, text="Preset Electrode").grid(row=1, column=0)
    probe_combo = ttk.Combobox(window, values=list(PRESETS), state="readonly")
    probe_combo.current(0)                       # default to "Custom (Manual)"
    probe_combo.grid(row=1, column=1)
    probe_combo.bind("<<ComboboxSelected>>", apply_preset)

    # Geometry fields, all in micrometres
    entries = []
    for i, (label, _) in enumerate(FIELDS):
        tk.Label(window, text=f"{label} (um)").grid(row=i + 2, column=0)
        entry = tk.Entry(window)
        entry.grid(row=i + 2, column=1)
        entries.append(entry)

    # Insertion site
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
