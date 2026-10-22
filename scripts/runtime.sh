#!/usr/bin/env bash
set -euo pipefail
project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
local_venv="${SHAV_VENV:-$HOME/.cache/shamonda-virus-phylogenomics/venv}"
export PATH="$local_venv/bin:$project_root/.venv/bin:$project_root/.runtime/bin:$PATH"
# Reuse scientific binaries without modifying the shared environment.
existing_runtime="/opt/homebrew/Caskroom/miniforge/base/envs/hantavirus-nextstrain/bin"
if [[ -d "$existing_runtime" ]]; then
  export PATH="$PATH:$existing_runtime"
fi
export MPLCONFIGDIR="$HOME/.cache/shamonda-virus-phylogenomics/matplotlib"
export PYTHONPYCACHEPREFIX="$HOME/.cache/shamonda-virus-phylogenomics/pycache"
export OMP_NUM_THREADS=2
cd -- "$project_root"
exec "$@"
