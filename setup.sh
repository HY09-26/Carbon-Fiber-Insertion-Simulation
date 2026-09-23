#!/usr/bin/env bash
#
# One-command setup:
#
#     bash setup.sh          # conda if available, otherwise pip
#     bash setup.sh --pip    # force pip + venv
#
# Builds the environment in ./.venv, registers the "Python (cf_v2)" Jupyter
# kernel the notebooks use, and checks that everything imports. Safe to re-run.

set -euo pipefail
cd "$(dirname "$0")"

USE_PIP=0
[ "${1:-}" = "--pip" ] && USE_PIP=1

ENV_DIR="./.venv"
PY="$ENV_DIR/bin/python"

# 1. Environment ---------------------------------------------------------------
if [ -x "$PY" ]; then
    echo "==> $ENV_DIR already exists, reusing it"
elif [ "$USE_PIP" -eq 0 ] && command -v conda >/dev/null 2>&1; then
    echo "==> Creating conda environment in $ENV_DIR"
    conda env create -p "$ENV_DIR" -f environment.yml
else
    echo "==> Creating venv in $ENV_DIR"
    python3 -m venv "$ENV_DIR"
    "$PY" -m pip install --quiet --upgrade pip
    "$PY" -m pip install --quiet -r requirements.txt
fi

# 2. Jupyter kernel named in the notebooks' metadata ---------------------------
echo "==> Registering Jupyter kernel 'cf_v2'"
"$PY" -m ipykernel install --user --name cf_v2 --display-name "Python (cf_v2)" >/dev/null

# 3. Verify --------------------------------------------------------------------
echo "==> Verifying"
"$PY" - <<'PYEOF'
import importlib, sys
missing = []
for m in ["numpy", "scipy", "pandas", "tifffile", "PIL", "pyvista", "ipykernel",
          "config", "Volume_bleeding", "Number_bleeding", "model_3D_visualization"]:
    try:
        importlib.import_module(m)
    except Exception as e:
        missing.append(f"{m} ({type(e).__name__}: {e})")
if missing:
    sys.exit("Setup incomplete:\n  " + "\n  ".join(missing))

import os, config
have = [a for a in config.ANIMALS if os.path.exists(config.tiff_path(a))]
print(f"    all imports OK; mouse_data has {len(have)}/{len(config.ANIMALS)} animals"
      + ("" if have else " - the tif stacks are not in the repo, ask the authors"))
PYEOF

cat <<'EOF'

Done. To reproduce the results:

  ./.venv/bin/ipython All_Simulation.ipynb         # electrode comparison, ~1 h
  ./.venv/bin/ipython Bleeding_per_diameter.ipynb  # cylinder control, several h

or open either notebook and choose the "Python (cf_v2)" kernel.
EOF
