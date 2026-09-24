# Carbon Fiber Insertion Simulation

Simulates inserting four neural probes into 3-D mouse cortical vasculature and
measures how much of the vascular network each one destroys, to test whether
thin carbon-fibre electrodes cause less vascular damage than commercial arrays.

| Probe | Key | Shank length / diameter (um) | Tip length / diameter (um) | Total (um) |
|---|---|---|---|---|
| Carbon Fiber 10 um | `CF` | 840 / 9.6 | 160 / 6.8 -> 0 | 1000 |
| Microprobes FMA | `FMA` | 971.91 / 25 | 28.09 / 25 -> 0 | 1000 |
| Shuttle 25 um | `Shuttle` | 990 / 25 | 10 / 25 (flat) | 1000 |
| Blackrock UEA | `UEA` | 950 / 90 -> 28 | 50 / 28 -> 3 | 1000 |

Each probe is a shank (near the surface) followed by a tip (at depth), both
tapering linearly, inserted straight down at 100 sites in each of 8 animals
(Reslice 2 and 3 are excluded: noisy raw data). Two metrics per insertion:

- **Volume** - vessel voxels inside the probe body.
- **Number** - distinct vessel components (3-D connected components) the probe touches.

A control sweep of plain cylinders, 4-80 um, separates the effect of thickness
from the effect of shape.

## Reproduce

```bash
git clone https://github.com/HY09-26/Carbon-Fiber-Insertion-Simulation.git
cd Carbon-Fiber-Insertion-Simulation
bash setup.sh
```

`setup.sh` builds the environment in `./.venv` (conda if available, pip
otherwise; ~2.3 GB) and registers the `Python (cf_v2)` Jupyter kernel.

Put the tif stacks in `mouse_data/` (~14 GB, not in the repo - ask the
authors), then run:

```bash
./.venv/bin/ipython All_Simulation.ipynb         # electrode comparison, ~1 h
./.venv/bin/ipython Bleeding_per_diameter.ipynb  # cylinder control, several hours
```

This runs a notebook headlessly and prints progress as it goes. Alternatively,
open it in Jupyter or VS Code and choose the `Python (cf_v2)` kernel.

About 16 GB of RAM is needed. The first run also segments each volume and
caches it in `seg_cache/`.

On Windows, run the steps in `setup.sh` by hand, using `.venv\Scripts\`
instead of `.venv/bin/`.

## Repository

```
config.py                    all settings: animals, probes, insertion protocol
All_Simulation.ipynb         electrode comparison -> results/electrodes/, results/figures/
Bleeding_per_diameter.ipynb  cylinder control     -> results/diameters/
Volume_bleeding.py           volume metric, probe geometry, data loading
Number_bleeding.py           count metric, vessel segmentation
model_3D_visualization.py    3-D rendering
UI.py                        GUI for exploring a single insertion (python UI.py)
setup.sh, environment.yml, requirements.txt
```

To add a probe or change the animals, edit `config.py`; both notebooks and the
GUI read from it. Lengths are in micrometres, which equal voxels in this data.

## Outputs

| File | Rows | Columns |
|---|---|---|
| `results/electrodes/Volume_<Key>.csv` | 800 (8 animals x 100 sites) | `file, position, overlap_area` |
| `results/electrodes/Number_<Key>.csv` | 800 | `file, position, overlap_number` |
| `results/diameters/<animal>_Volume.csv` | 100 sites | `position, 4um, 8um, ..., 80um` |
| `results/diameters/<animal>_Number.csv` | 100 sites | `position, 4um, 8um, ..., 80um` |
| `results/figures/<Key>_<animal>.tif` | - | 3-D render at one site, illustrative (not tracked) |
