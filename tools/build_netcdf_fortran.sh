#!/usr/bin/env bash
# -----------------------------------------------------------------------
# build_netcdf_fortran.sh — build netcdf-fortran for ONE Fortran compiler
#
# Rakali's only compiler-coupled dependency is netcdf-fortran. Fortran
# `.mod` files and module symbol manglings are compiler-specific, so a
# netcdf-fortran built by gfortran cannot be consumed by nvfortran, ifx or
# flang -- you need one per compiler.
#
# Everything BELOW it -- netcdf-c, HDF5, zlib -- is C, and a C library is
# reusable across all of these compilers on Linux x86-64. So you do NOT
# need a per-compiler netcdf-c: take a prebuilt one from your distro,
# module tree, or conda, and build only this Fortran wrapper on top. That
# is minutes of work instead of an HDF5-from-source stack.
#
# No netcdf-c anywhere (only HDF5)? Use tools/build_netcdf.sh, which
# builds netcdf-c too. This script is that one's --fortran-only mode.
#
# Usage:
#   tools/build_netcdf_fortran.sh --fc nvfortran
#   tools/build_netcdf_fortran.sh --fc ifx --prefix ~/opt/netcdff-ifx
#   tools/build_netcdf_fortran.sh --fc gfortran --netcdf-c /usr
#
# netcdf-c is located in this order:
#   1. --netcdf-c <prefix>
#   2. $NETCDFC_DIR
#   3. `nc-config --prefix`, if nc-config is on PATH
#
# Other options: --prefix DIR, --version V (netcdf-fortran release),
# --jobs N, --keep-build; anything else is passed to build_netcdf.sh
# (e.g. --offline DIR, --dry-run, --static, --cc CC).
#
# Note that Rakali needs NEITHER parallel netcdf NOR HDF5's Fortran
# bindings: its I/O is per-rank serial `nf90_create` with an offline merge
# (tools/merge_output.py), and it imports `netcdf` and nothing else. A
# serial netcdf-c against a serial `hdf5 ~fortran` is enough.
# -----------------------------------------------------------------------
set -euo pipefail

args=()
while [[ $# -gt 0 ]]; do
   case "$1" in
      --version) args+=(--netcdf-fortran-version "${2:-}"); shift 2 ;;
      -h|--help) sed -n '3,/^# ------/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit 0 ;;
      *)         args+=("$1"); shift ;;
   esac
done

exec "$(dirname "$0")/build_netcdf.sh" --fortran-only "${args[@]}"
