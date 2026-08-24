#!/usr/bin/env bash
#
# One-command setup for CF_v2.
#
#     bash setup.sh              # conda if available, otherwise pip
#     bash setup.sh --pip        # force pip + venv
#     bash setup.sh --exact      # pin every version (requirements-lock.txt)
#
# Creates the environment in ./.venv, registers a "Python (cf_v2)" Jupyter
# kernel, writes .vscode/settings.json, and verifies that everything imports.
#
# Windows: see README.md, "Setting up the environment".

set -euo pipefail
cd "$(dirname "$0")"

USE_PIP=0
EXACT=0
for arg in "$@"; do
    case "$arg" in
        --pip)   USE_PIP=1 ;;
        --exact) EXACT=1 ;;
        -h|--help) sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown option: $arg (try --help)" >&2; exit 1 ;;
    esac
done

ENV_DIR="./.venv"
PY="$ENV_DIR/bin/python"

# ---------------------------------------------------------------------------
# 1. Create the environment
# ---------------------------------------------------------------------------
if [ -x "$PY" ]; then
    echo "==> $ENV_DIR already exists, reusing it"
elif [ "$USE_PIP" -eq 0 ] && command -v conda >/dev/null 2>&1; then
    echo "==> Creating conda environment in $ENV_DIR"
    # -p puts the environment inside the project, so the VS Code settings
    # below work unchanged on every machine.
    conda env create -p "$ENV_DIR" -f environment.yml
else
    if [ "$USE_PIP" -eq 0 ]; then
        echo "==> conda not found, falling back to pip"
    fi
    echo "==> Creating venv in $ENV_DIR"
    python3 -m venv "$ENV_DIR"
    "$PY" -m pip install --quiet --upgrade pip
fi

# ---------------------------------------------------------------------------
# 2. Install / pin dependencies
# ---------------------------------------------------------------------------
if [ "$EXACT" -eq 1 ]; then
    echo "==> Installing exact pinned versions (requirements-lock.txt)"
    "$PY" -m pip install --quiet -r requirements-lock.txt
elif [ ! -f "$ENV_DIR/conda-meta/history" ]; then
    # pip route: conda env create already installed everything.
    echo "==> Installing dependencies (requirements.txt)"
    "$PY" -m pip install --quiet -r requirements.txt
fi

# ---------------------------------------------------------------------------
# 3. Register the Jupyter kernel the notebooks ask for
# ---------------------------------------------------------------------------
echo "==> Registering Jupyter kernel 'cf_v2'"
"$PY" -m ipykernel install --user --name cf_v2 --display-name "Python (cf_v2)" >/dev/null

# ---------------------------------------------------------------------------
# 4. VS Code settings (never overwrite an existing one)
# ---------------------------------------------------------------------------
mkdir -p .vscode
if [ -f .vscode/settings.json ]; then
    echo "==> .vscode/settings.json already exists, leaving it alone"
else
    cp .vscode/settings.json.example .vscode/settings.json
    echo "==> Wrote .vscode/settings.json"
fi

# ---------------------------------------------------------------------------
# 5. Verify
# ---------------------------------------------------------------------------
echo "==> Verifying"
"$PY" - <<'PYEOF'
import importlib, sys, os
print(f"    python     {sys.version.split()[0]}  ({sys.executable})")
missing = []
for m in ["numpy", "scipy", "pandas", "tifffile", "PIL", "pyvista",
          "tqdm", "matplotlib", "ipykernel"]:
    try:
        mod = importlib.import_module(m)
        print(f"    {m:11s} {getattr(mod, '__version__', 'ok')}")
    except ImportError:
        missing.append(m)
        print(f"    {m:11s} MISSING")
sys.path.insert(0, ".")
for m in ["Volume_bleeding", "Number_bleeding", "model_3D_visualization"]:
    try:
        importlib.import_module(m)
    except Exception as e:
        missing.append(m)
        print(f"    {m}: {type(e).__name__}: {e}")
if missing:
    sys.exit(f"\nSetup incomplete, missing: {', '.join(missing)}")
n = len(os.listdir("mouse_data")) if os.path.isdir("mouse_data") else 0
print(f"\n    mouse_data/: {n} files"
      + ("" if n else "  <- empty; see README.md, 'Getting the data'"))
PYEOF

cat <<'EOF'

Done.

  Notebook   open All_Simulation.ipynb and pick the kernel "Python (cf_v2)"
  Shell      ./.venv/bin/python UI.py
  Data       put the tif stacks in mouse_data/ (see README.md)
EOF
