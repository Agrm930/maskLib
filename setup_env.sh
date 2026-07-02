#!/usr/bin/env bash
# =============================================================================
#  maskLib Environment Setup Script
#  Sets up a conda environment and installs all dependencies for maskLib.
#  Run this script from the root of the maskLib repository.
# =============================================================================

set -e

ENV_NAME="masklib"
PYTHON_VERSION="3.10"

echo "=============================================="
echo "  maskLib Environment Setup"
echo "=============================================="

# ---- 1. Create conda environment ----
echo ""
echo "[1/4] Creating conda environment: $ENV_NAME (Python $PYTHON_VERSION)"
conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y

# Activate
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# ---- 2. Install core Python dependencies ----
echo ""
echo "[2/4] Installing Python dependencies..."

# Core dependencies (always required)
pip install dxfwrite          # DXF file writer — the main drawing engine
pip install numpy             # Numerical arrays
pip install matplotlib        # Plotting / path utilities (used in waffle functions)

# Optional but commonly used
pip install ezdxf             # Modern DXF reader/writer (used in some modules)
pip install gdspy             # GDSII layout (optional, used in some modules)
pip install opencv-python     # cv2 (optional, used in some modules)

# ---- 3. Install maskLib in editable mode ----
echo ""
echo "[3/4] Installing maskLib in editable mode..."
pip install -e .

echo ""
echo "[4/4] Verifying installation..."
python -c "
import maskLib
from maskLib import MaskLib, microwaveLib, fluxoniumLib, junctionLib
from maskLib import qubitLib, resonatorLib, markerLib, utilities
print('  All maskLib modules imported successfully.')
print('  maskLib is ready to use!')
"

echo ""
echo "=============================================="
echo "  Setup complete!"
echo "  Activate your environment with:"
echo "    conda activate $ENV_NAME"
echo "=============================================="
