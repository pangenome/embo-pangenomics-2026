#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$root_dir/.cache/xdg" "$root_dir/.cache/quarto"

export XDG_CACHE_HOME="$root_dir/.cache/xdg"
export QUARTO_CACHE_DIR="$root_dir/.cache/quarto"

exec quarto render "$@"
