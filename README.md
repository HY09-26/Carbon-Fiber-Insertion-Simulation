# CF_v2 - Electrode Insertion / Vascular Damage Simulation

Simulates inserting five different neural probes into 3-D mouse cortical
vasculature and measures how much of the vascular network each one destroys.
The motivating question is whether thin, tapered carbon-fibre electrodes cause
less vascular damage than the larger commercial arrays.

Probes compared:

| Probe | Shank length / diameter (um) | Tip length / diameter (um) | Total (um) |
|---|---|---|---|
| Carbon Fiber | 840 / 8.4 | 160 / 6.8 -> 0 | 1000 |
| Paradromics Connexus | 875 / 20 | 125 / 20 -> 0 | 1000 |
| Microprobes FMA | 971.91 / 25 | 28.09 / 25 -> 0 | 1000 |
| Neuralink Shuttle | 990 / 25 | 10 / 25 (flat) | 1000 |
| Blackrock UEA | 950 / 90 -> 28 | 50 / 28 -> 3 | 1000 |

Each probe is modelled as two stacked conical frusta: a **shank** near the
cortical surface and a **tip** at depth, both tapering linearly. The insertion
path is a straight line perpendicular to the cortical surface.

Two damage metrics are computed:

- **Volume** - number of vessel *voxels* swept by the probe body.
- **Number** - number of distinct *connected components* of the segmented
  vasculature the probe touches.

A control sweep of plain cylinders (diameter 4-80 um, step 4) separates the
effect of raw thickness from the effect of taper and shape.

---

## Quick start

```bash
bash setup.sh
```

That is the whole setup. It creates the environment in `./.venv`, installs
everything, registers a `Python (cf_v2)` Jupyter kernel and writes
`.vscode/settings.json` - no paths to edit. conda is used if it is on your
PATH, otherwise it falls back to pip.

Then put the raw data in `mouse_data/` (see [Getting the data](#getting-the-data)),
open `All_Simulation.ipynb`, pick the **Python (cf_v2)** kernel, and run it top
to bottom.

---

## Requirements


### Software

Python 3.10 or newer. Dependencies are declared three ways, pick one:

| File | Use it when |
|---|---|
| `environment.yml` | **Recommended.** conda-forge ships prebuilt VTK binaries for macOS (arm64 + x86_64), Linux and Windows, so pyvista installs cleanly. |
| `requirements.txt` | No conda available. Version floors only, so it stays portable. |
| `requirements-lock.txt` | You want to reproduce the published numbers against the exact versions they were computed with. |

`tkinter` (used only by `UI.py`) ships with CPython. On some Linux
distributions it is a separate system package, e.g. `apt install python3-tk`.

---

## Setting up the environment

`bash setup.sh` does everything in this section. Read on only if you want to
do it by hand, or you are on Windows.

| Command | Effect |
|---|---|
| `bash setup.sh` | conda if it is on your PATH, otherwise pip |
| `bash setup.sh --pip` | force pip + venv even when conda is installed |
| `bash setup.sh --exact` | pin every version from `requirements-lock.txt` |

Re-running is safe: an existing `./.venv` is reused, and an existing
`.vscode/settings.json` is never overwritten.

The environment **must** live at `./.venv`, because
`.vscode/settings.json.example` points `python.defaultInterpreterPath` at
`${workspaceFolder}/.venv/bin/python`. That is what makes the VS Code settings
work on any machine without editing.

### By hand (macOS / Linux)

```bash
# conda route
conda env create -p ./.venv -f environment.yml
# or pip route
python3 -m venv .venv && ./.venv/bin/pip install -r requirements.txt

./.venv/bin/python -m ipykernel install --user --name cf_v2 --display-name "Python (cf_v2)"
cp .vscode/settings.json.example .vscode/settings.json
```

### Windows

```powershell
conda env create -p .\.venv -f environment.yml   # or: python -m venv .venv
.\.venv\Scripts\pip install -r requirements.txt  # pip route only
.\.venv\Scripts\python -m ipykernel install --user --name cf_v2 --display-name "Python (cf_v2)"
copy .vscode\settings.json.example .vscode\settings.json
```

Then edit `.vscode/settings.json` and change `bin/python` to
`Scripts/python.exe`.

### Notes

`.vscode/settings.json` is gitignored; `.vscode/settings.json.example` is the
version under version control.

`jupyter.notebookFileRoot` is set to `${fileDirname}` because cell 0 of each
notebook resolves `mouse_data` relative to the working directory, so a
notebook must run from this folder.

If your machine has other Python installations that advertise themselves as
kernels, list them under `jupyter.kernels.filter` to keep them out of the
kernel picker. Entries that do not exist on your machine are ignored.

---


## Running the pipeline

### 1. `All_Simulation.ipynb` - main results

| Cell | What it does |
|---|---|
| 0 | Imports; builds `tiff_files` from `mouse_data/` (8 animals) |
| 1 | All parameters: insertion site, crop margins, `pos_num = 100`, `depth_limit = 1000`, and the five probe geometries |
| 2 | Segments every volume into connected components and caches the result in `seg_cache/*.npz`. Skipped automatically if the cache exists |
| 3-7 | One cell per probe (UEA, Neuralink, FMA, Paradromics, Carbon Fiber). Each runs both metrics over all 8 animals, writes two CSVs, and renders one 3-D figure per animal |

Cells 3-7 are five copies of the same code differing only in
`selected_electrode`. Folding them into a loop is the obvious cleanup but has
been left alone to avoid changing behaviour.

**Runtime**, measured on an Apple-silicon laptop:

- loading one volume: ~2 s (3.0 GB resident)
- loading one cached segmentation: ~5 s (6.0 GB resident)
- ~0.5 s per insertion per metric
- one probe across 8 animals: **~13 min**; all five probes: **~1 hour**, plus
  3-D rendering
- the first run also pays for segmentation, which is not cached yet

### 2. `Bleeding_per_diameter.ipynb` - cylinder control

Sweeps plain cylinders of diameter 4, 8, ..., 80 um (`np.arange(4, 81, 4)`)
over the same 8 animals and 100 sites, giving a dose-response curve against
thickness alone. Reuses `seg_cache/`.

Expect **several hours**: 20 diameters x 100 sites x 8 animals x 2 metrics.

### 3. `UI.py` - single insertion, interactive

```bash
python UI.py
```

Pick a tif, choose a probe preset (or type custom dimensions), set an
insertion site, and get the vessel-voxel count plus the 3-D view. Exploration
tool only - the paper's numbers come from the notebooks.

Note that the GUI takes **radii**, while the underlying functions take
**diameters**; the conversion happens inside `run_simulation`.

---

## Outputs

```
output_csv/
  Volume_<Probe>.csv    columns: file, position, overlap_area
  Number_<Probe>.csv    columns: file, position, overlap_number
                        1 row per (animal, insertion site) = 800 rows

Volume_bleeding_per_diameter/
  <animal>_Volume.csv   columns: position, 4um, 8um, ... 80um
Number_bleeding_per_diameter/
  <animal>_Number.csv   columns: position, 4um, 8um, ... 80um
                        1 row per insertion site = 100 rows

output_images/
  <Probe>_<animal>.tif  3-D render, 300 dpi
                        gold = vessels outside the probe
                        red  = vessels inside the probe
                        grey = the probe

seg_cache/
  <animal>_seg.npz      cached connected-component labelling (int32)
```

`seg_cache/` is derived data - delete it and cell 2 regenerates it. It has
been checked against the current raw tifs: the only difference is 2,621
single-voxel components dropped by the `min_size=2` filter.

`seg_cache/` and `output_images/` are gitignored; the result CSVs in
`output_csv/` and the two `*_bleeding_per_diameter/` folders are tracked.

---

## Code guide

```
Volume_bleeding.py           volume metric  (core)
Number_bleeding.py           count metric + segmentation (core)
model_3D_visualization.py    3-D rendering
Area_UI.py                   deprecated shim -> Volume_bleeding
UI.py                        Tkinter GUI for a single insertion
All_Simulation.ipynb         main pipeline, 5 probes x 8 animals
Bleeding_per_diameter.ipynb  cylinder control sweep
```

### `Volume_bleeding.py`

The volume metric, and the geometry primitives everything else reuses.

| Function | Purpose |
|---|---|
| `create_cone_mask(height, width, x_center, y_center, radius)` | Boolean disc for one depth slice. Uses `np.ogrid`, so no full 2-D coordinate grid is built. |
| `calculate_cone_radius(depth, base_diameter, top_diameter, total_length)` | Radius of a linearly tapering segment. `base == top` gives a cylinder. |
| `simulate_cone_insertion(...)` | One insertion. Walks down from `start_slice`, picks shank or tip geometry per depth, sums vessel voxels inside the disc. Stops once past `shank_length + tip_length`. |
| `find_first_black_pixel_slice(img_data, x_center, y_center)` | Cortical surface depth: first voxel equal to 0 down a column. |
| `process_cone_positions(...)` | Repeats the insertion at `pos_num` evenly spaced x positions, `space = height / (pos_num + 1)` so nothing lands on the edge. |

### `Number_bleeding.py`

The count metric, plus the segmentation that makes it possible.

| Function | Purpose |
|---|---|
| `label_vessel(img, min_size, connectivity)` | `scipy.ndimage.label` into 3-D connected components, then drops components smaller than `min_size`. `connectivity=2` means a 3x3x3 structuring element (26-connectivity). |
| `vessel_seg(img_data, ...)` | Thin wrapper used by the notebooks. The `merge_labels` step is commented out. |
| `merge_labels(labeled_data, distance)` | **Not used.** Would merge labels whose surfaces come within `distance` voxels, to repair vessels broken by imperfect segmentation. Bounding-box rejection, then a KD-tree on surface voxels, then union-find. The pairwise loop is O(n^2) over ~41,000 labels, so it does not finish. |
| `get_all_surface_voxels_vectorized(...)` | Surface voxels per label, via one binary erosion. Helper for `merge_labels`. |
| `simulate_cone_insertion_num(...)` | Same geometry walk as the volume version, but collects label ids into a set, so a vessel crossed at several depths counts once. Has the bounds guard the volume version lacks. |
| `process_cone_positions_num(...)` | The `pos_num`-site loop for the count metric. Takes the **label image**, not the raw binary volume. |
| `get_top_n_labels`, `visualize_labeled_3d` | Debug helpers for eyeballing the segmentation. |

Note that this module also imports `calculate_cone_radius` and
`find_first_black_pixel_slice` from elsewhere and then re-defines both. The
local definitions win. The only difference is that the local
`find_first_black_pixel_slice` returns 0 instead of `None` for a column that
never reaches background.

### `model_3D_visualization.py`

| Function | Purpose |
|---|---|
| `cropping_img(...)` | Cuts a slab around the insertion (default 160 x 100 x 1000) so the renderer stays responsive, and re-expresses the insertion centre in cropped coordinates. |
| `visualize_cone_pyvista(...)` | Builds a scalar field (1 = vessel outside probe, 2 = vessel inside), thresholds each into a mesh, adds a 36-gon wireframe cone for the probe, a bounding box and a 100 um scale bar. `ui=1` opens an interactive window; `ui=0` returns the plotter for a headless `.screenshot()`. Key bindings: `s` screenshot, `r` reset camera. |

### `Area_UI.py`

Deprecated. It used to hold a second copy of the volume-metric code, identical
to `Volume_bleeding.py` apart from one progress `print`. It now re-exports
from `Volume_bleeding` so the algorithm lives in one place. Import
`Volume_bleeding` directly in new code.

### `UI.py`

Tkinter GUI described under [Running the pipeline](#3-uipy--single-insertion-interactive).
It inlines its own copy of the crop logic instead of calling `cropping_img`;
there is a `TODO` on that line.

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'pyvista'`**
The notebook is running on the wrong kernel. Pick the kernel that belongs to
the environment you created (top right of the notebook in VS Code, or
Kernel > Change kernel in Jupyter). If several interpreters on your machine
look alike, use `jupyter.kernels.filter` in `.vscode/settings.json` to hide
the wrong ones.

**`AssertionError: ... not found - run this notebook from the CF_v2 folder`**
Cell 0 resolves `mouse_data` relative to the working directory. Set
`"jupyter.notebookFileRoot": "${fileDirname}"` in VS Code, or start Jupyter
from this folder.

**Rendering fails on a headless machine**
`pyvista` needs a display. Either start a virtual framebuffer before
rendering:

```python
import pyvista as pv
pv.start_xvfb()
```

or set `PYVISTA_OFF_SCREEN=true` in the environment. Both need `xvfb`
installed.

**`MemoryError` or heavy swapping during cell 2**
The segmentation pass is the peak: the `uint16` volume is 3.0 GB and the
`int32` label image is 6.0 GB, held together with SciPy's internals. 16 GB is
the practical floor. If the cache in `seg_cache/` is already populated, cell 2
only loads it and the peak is lower.

**`numpy` fails to import with an architecture error on macOS**
An x86_64 numpy in an arm64 interpreter (or the reverse). Build the
environment from `environment.yml` rather than reusing a system Python.

---

## Conventions

- All lengths, diameters and radii are in **micrometres**, which equals voxels
  at this 1 um isotropic sampling.
- The simulation functions take **diameters**; the GUI takes **radii** and
  doubles them.
- Volumes must be transposed to `(z, x, y)` with
  `np.transpose(img, axes=(0, 2, 1))` before being passed to any function here.
- The whole repository is ASCII-only. The micrometre unit is written `um`, not
  the Greek mu.
