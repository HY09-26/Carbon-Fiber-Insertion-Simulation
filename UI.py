"""
UI.py
=====
Small Tkinter front-end for a SINGLE insertion.

Pick a .tif volume, choose one of the five modelled probes (or type your own
dimensions), and run one insertion at a chosen (x, y).  It reports the
intersected vessel-voxel count and opens the 3-D view.

This is an exploration tool.  The numbers that go into the paper come from
`All_Simulation.ipynb`, which sweeps 100 insertion sites across every animal.

Run it with:

    python UI.py

Dimensions are entered as RADII in micrometres; they are doubled to diameters
before being handed to the simulation functions, which work in diameters.
"""

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import tifffile as tiff
import numpy as np


# Importing functions from another file
# NOTE: the historical version also imported the same two names from
# `Area_UI`, which the `Volume_bleeding` import below then shadowed.
# `Area_UI` is now a deprecated re-export of `Volume_bleeding`, so that
# redundant line has been removed.
from model_3D_visualization import visualize_cone_pyvista, calculate_cone_radius, cropping_img
from Volume_bleeding import simulate_cone_insertion, find_first_black_pixel_slice, process_cone_positions

## Default Electrodes
# Manufacturer dimensions, in micrometres. "tip" is the deepest segment,
# "shank" the part nearer the cortical surface; both taper linearly from
# radius_start to radius_end. Shank + tip length is 1000 um for every probe.
preset_configs = {
    "Custom (Manual)": {},
    "Carbon Fiber (8.4um)": {
        "tip_length": 160, "tip_radius_start": 3.4, "tip_radius_end": 0,
        "shank_length": 840, "shank_radius_start": 4.2, "shank_radius_end": 4.2
    },
    "Paradromics Connexus (20um)": {
        "tip_length": 125, "tip_radius_start": 10, "tip_radius_end": 0,
        "shank_length": 875, "shank_radius_start": 10, "shank_radius_end": 10
    },
    "Microprobes FMA (25um)": {
        "tip_length": 28.09, "tip_radius_start": 12.5, "tip_radius_end": 0,
        "shank_length": 971.91, "shank_radius_start": 12.5, "shank_radius_end": 12.5
    },
    "Neuralink Shuttle (25um)": {
        "tip_length": 10, "tip_radius_start": 12.5, "tip_radius_end": 12.5,
        "shank_length": 990, "shank_radius_start": 12.5, "shank_radius_end": 12.5
    },
    "Blackrock UEA (90um)": {
        "tip_length": 50, "tip_radius_start": 14, "tip_radius_end": 1.5,
        "shank_length": 950, "shank_radius_start": 45, "shank_radius_end": 14
    }
}


def run_gui():
    """Build and start the Tkinter window. Blocks until the window is closed."""
    print("UI launching...")

    def browse_file():
        """File picker for the .tif volume."""
        file_path = filedialog.askopenfilename(filetypes=[("TIFF files", "*.tif *.tiff")])
        entry_file.delete(0, tk.END)
        entry_file.insert(0, file_path)

    def apply_preset(event=None):
        """Fill the dimension fields from the selected probe preset."""
        selection = electrode_combo.get()
        config = preset_configs.get(selection, {})
        if config:
            entry_tip_length.delete(0, tk.END)
            entry_tip_length.insert(0, str(config["tip_length"]))
            entry_tip_radius_start.delete(0, tk.END)
            entry_tip_radius_start.insert(0, str(config["tip_radius_start"]))
            entry_tip_radius_end.delete(0, tk.END)
            entry_tip_radius_end.insert(0, str(config["tip_radius_end"]))
            entry_shank_length.delete(0, tk.END)
            entry_shank_length.insert(0, str(config["shank_length"]))
            entry_shank_radius_start.delete(0, tk.END)
            entry_shank_radius_start.insert(0, str(config["shank_radius_start"]))
            entry_shank_radius_end.delete(0, tk.END)
            entry_shank_radius_end.insert(0, str(config["shank_radius_end"]))

    def run_simulation():
        """Load the volume, run one insertion, report the count, show the 3-D view."""
        try:
            filepath = entry_file.get()
            img_data = tiff.imread(filepath)
            # (z, y, x) as stored on disk -> (z, x, y) as the simulation expects.
            img_data = np.transpose(img_data, axes=(0, 2, 1))
            img_data = img_data.astype(np.uint16)

            print(img_data.shape)

            tip_length = float(entry_tip_length.get())
            tip_radius_start = float(entry_tip_radius_start.get())
            tip_radius_end = float(entry_tip_radius_end.get())
            shank_length = float(entry_shank_length.get())
            shank_radius_start = float(entry_shank_radius_start.get())
            shank_radius_end = float(entry_shank_radius_end.get())
            x_center = int(entry_x.get())
            y_center = int(entry_y.get())

            # === Find the starting slice ===
            start_slice = find_first_black_pixel_slice(img_data, x_center, y_center)

            # Radii are doubled because the simulation takes diameters.
            result = simulate_cone_insertion(
                img_data, x_center, y_center,
                shank_length, shank_radius_start * 2, shank_radius_end * 2,
                tip_length, tip_radius_start * 2, tip_radius_end * 2,
                start_slice, 1000
            )


            # === Cropping ===
            # Same slab as model_3D_visualization.cropping_img, inlined here.
            # TODO: call cropping_img() instead, so the crop logic lives in one place.
            crop_margin_x = 80
            crop_margin_y = 50
            x_min = max(x_center - crop_margin_x , 0)
            x_max = min(x_center + crop_margin_x , img_data.shape[1])
            y_min = y_center - crop_margin_y
            y_max = y_center + crop_margin_y
            cropped_img_data = img_data[start_slice:start_slice+1000, x_min:x_max, y_min:y_max]  # shape: [depth, H, W]


            # ===  the cone's center to match the cropped region ===
            adjusted_x_center = crop_margin_x
            adjusted_y_center = crop_margin_y


            # === Visualize using the cropped data ===
            visualize_cone_pyvista(
                cropped_img_data, adjusted_x_center, adjusted_y_center, shank_length=shank_length,
                shank_base_diameter=shank_radius_start * 2, shank_top_diameter=shank_radius_end * 2,
                tip_length=tip_length, tip_base_diameter=tip_radius_start * 2, tip_top_diameter=tip_radius_end * 2,
                start_slice=start_slice, depth_limit=1000, ui=1
            )

            messagebox.showinfo("Result", f"Total vessel voxels intersected: {int(result)}")
        except Exception as e:
            messagebox.showerror("Error", str(e))




    window = tk.Tk()
    window.title("Visualization:Electrode Insertion Simulation ")

    # File input
    tk.Label(window, text="Load File").grid(row=0, column=0)
    entry_file = tk.Entry(window, width=40)
    entry_file.grid(row=0, column=1)
    tk.Button(window, text="Browse", command=browse_file).grid(row=0, column=2)

    # Electrode Preset Dropdown
    tk.Label(window, text="Preset Electrode").grid(row=1, column=0)
    electrode_combo = ttk.Combobox(window, values=list(preset_configs.keys()), state="readonly")
    electrode_combo.current(0)  # Default to "Custom"
    electrode_combo.grid(row=1, column=1)
    electrode_combo.bind("<<ComboboxSelected>>", apply_preset)

    # Parameter input (all lengths and radii in micrometres)
    labels = [
        "Tip Length", "Tip Radius Start", "Tip Radius End",
        "Shank Length", "Shank Radius Start", "Shank Radius End",
        "X Center", "Y Center"
    ]
    entries = []
    for i, label in enumerate(labels):
        tk.Label(window, text=label).grid(row=i+2, column=0)
        entry = tk.Entry(window)
        entry.grid(row=i+2, column=1)
        entries.append(entry)

    # Unpack in the same order as `labels` above.
    (
        entry_tip_length, entry_tip_radius_start, entry_tip_radius_end,
        entry_shank_length, entry_shank_radius_start, entry_shank_radius_end,
        entry_x, entry_y
    ) = entries

    # Simulation button
    tk.Button(window, text="Run Simulation", command=run_simulation).grid(row=len(labels)+3, column=1)

    window.mainloop()

if __name__ == "__main__":
    run_gui()
