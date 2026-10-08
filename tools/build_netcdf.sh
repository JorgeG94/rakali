#!/usr/bin/env bash
# -----------------------------------------------------------------------
# build_netcdf.sh -- build netcdf-c AND netcdf-fortran into one prefix
#
# For machines that have a compiler, (optionally) MPI and an HDF5, but no
# netcdf-c you can use -- e.g. NERSC Perlmutter (cray-hdf5 + PrgEnv-*).
# It downloads pinned, checksummed Unidata release tarballs, builds
# netcdf-c against YOUR HDF5 (it never builds HDF5), builds netcdf-fortran
# with YOUR Fortran compiler on top, smoke-tests the pair, writes
# <prefix>/rakali-netcdf.env and prints the cmake line for Rakali.
#
# Usage:
#   tools/build_netcdf.sh --fc gfortran --prefix ~/opt/netcdf-gnu
#   tools/build_netcdf.sh --cc cc --fc ftn --prefix $SCRATCH/netcdf-gnu  # Cray
#   tools/build_netcdf.sh --cc nvc --fc nvfortran --hdf5 /opt/hdf5 --prefix P
#
#   # no internet on the build node: fetch on a login node, build elsewhere
#   tools/build_netcdf.sh --download-only --download-dir ~/netcdf-src
#   tools/build_netcdf.sh --offline ~/netcdf-src --fc ftn --prefix P
#
# Options:
#   --fc FC               Fortran compiler (default $FC; ftn on a Cray PE)
#   --cc CC               C compiler for netcdf-c (default $CC, else the one
#                         paired with --fc: gfortran->gcc, nvfortran->nvc,
#                         ifx->icx, flang->clang, ftn->cc)
#   --hdf5 PREFIX         HDF5 install to build against. Else $HDF5_DIR,
#                         $HDF5_ROOT, then `h5cc -show` (`h5pcc` with --mpi),
#                         then pkg-config. HDF5 is never built here.
#   --mpi                 build against a PARALLEL HDF5 with CC=mpicc,
#                         FC=mpifort (cc/ftn on Cray). Rakali itself only
#                         needs SERIAL netcdf -- use this only when the HDF5
#                         you have is parallel. Default: serial.
#   --prefix DIR          install prefix (default ~/opt/netcdf-<ver>-<fc>)
#   --jobs N, -j N        parallel make jobs (default: nproc)
#   --netcdf-c-version V  default 4.9.3 (checksum pinned)
#   --netcdf-fortran-version V
#                         default 4.6.2 (checksum pinned)
#   --netcdf-c-sha256 S / --netcdf-fortran-sha256 S
#                         checksum for a version not in the built-in table
#   --static              static libraries only (older Cray PEs link static)
#   --cflags "..." / --fflags "..."
#                         extra compiler flags, e.g. "-tp=zen3" for NVHPC
#   --check               also build + run the upstream test suites
#   --download-dir DIR    keep/reuse verified tarballs in DIR
#   --download-only       fetch + verify tarballs into --download-dir, stop
#   --offline DIR         use tarballs from DIR, never touch the network
#   --work-dir DIR        build tree location (default: mktemp under $TMPDIR)
#   --keep-build          keep the build tree (and its logs) afterwards
#   --dry-run             resolve everything, print the commands, build nothing
#   --fortran-only        skip netcdf-c; build netcdf-fortran on an EXISTING
#                         netcdf-c (--netcdf-c PREFIX, else $NETCDFC_DIR, else
#                         `nc-config --prefix`). tools/build_netcdf_fortran.sh
#                         is this mode.
#   -h, --help            this text
#
# netcdf-c is built with its autotools `configure`, not CMake: it takes the
# Cray wrappers (cc/ftn) and NVHPC (nvc) as plain $CC with no toolchain
# file, finds HDF5 from CPPFLAGS/LDFLAGS whatever its layout, and does not
# depend on the CMake FindHDF5 module guessing right. netcdf-fortran uses
# its CMake build (what CI has exercised for every compiler).
# -----------------------------------------------------------------------
set -euo pipefail

# --- pinned releases ----------------------------------------------------
NC_VERSION="4.9.3"
NCF_VERSION="4.6.2"
URL_BASE="https://downloads.unidata.ucar.edu"

# sha256 of the Unidata release tarballs (downloads.unidata.ucar.edu).
known_sha256() {
   case "$1" in
      netcdf-c-4.9.3)       echo a474149844e6144566673facf097fea253dc843c37bc0a7d3de047dc8adda5dd ;;
      netcdf-c-4.10.1)      echo db3b69ff4a5ee1a7d79a5c36664d2128b752c266e966369fcf7311ec5f927564 ;;
      netcdf-fortran-4.6.2) echo df26b99d9003c93a8bc287b58172bf1c279676f8c10d6dd0daf8bc7204877096 ;;
      netcdf-fortran-4.6.4) echo 98159c1e0f63b3b59bb5eda12f2d80126f5b1aad93032d1490989a5752e0df99 ;;
      *) echo "" ;;
   esac
}

# --- defaults -----------------------------------------------------------
FC="${FC:-}"
CC_ARG=""
HDF5_ARG=""
USE_MPI=0
PREFIX=""
JOBS="${JOBS:-$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
NC_SHA256=""
NCF_SHA256=""
STATIC=0
EXTRA_CFLAGS=""
EXTRA_FFLAGS=""
CHECK=0
DOWNLOAD_DIR=""
DOWNLOAD_ONLY=0
OFFLINE=0
WORK_ARG=""
KEEP_BUILD=0
DRY_RUN=0
FORTRAN_ONLY=0
NETCDF_C=""

die() { printf '\nerror: %s\n' "$*" >&2; exit 1; }
note() { printf '%s\n' "$*"; }
warn() { printf 'warning: %s\n' "$*" >&2; }

usage() {
   sed -n '3,/^# ------/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'
   exit "${1:-0}"
}

need_arg() { [[ $# -ge 2 && -n "$2" ]] || die "$1 needs a value"; }

while [[ $# -gt 0 ]]; do
   case "$1" in
      --fc)                     need_arg "$@"; FC="$2"; shift 2 ;;
      --cc)                     need_arg "$@"; CC_ARG="$2"; shift 2 ;;
      --hdf5)                   need_arg "$@"; HDF5_ARG="$2"; shift 2 ;;
      --mpi)                    USE_MPI=1; shift ;;
      --prefix)                 need_arg "$@"; PREFIX="$2"; shift 2 ;;
      --jobs|-j)                need_arg "$@"; JOBS="$2"; shift 2 ;;
      --netcdf-c-version)       need_arg "$@"; NC_VERSION="$2"; shift 2 ;;
      --netcdf-fortran-version) need_arg "$@"; NCF_VERSION="$2"; shift 2 ;;
      --netcdf-c-sha256)        need_arg "$@"; NC_SHA256="$2"; shift 2 ;;
      --netcdf-fortran-sha256)  need_arg "$@"; NCF_SHA256="$2"; shift 2 ;;
      --static)                 STATIC=1; shift ;;
      --cflags)                 need_arg "$@"; EXTRA_CFLAGS="$2"; shift 2 ;;
      --fflags)                 need_arg "$@"; EXTRA_FFLAGS="$2"; shift 2 ;;
      --check)                  CHECK=1; shift ;;
      --download-dir)           need_arg "$@"; DOWNLOAD_DIR="$2"; shift 2 ;;
      --download-only)          DOWNLOAD_ONLY=1; shift ;;
      --offline)                need_arg "$@"; OFFLINE=1; DOWNLOAD_DIR="$2"; shift 2 ;;
      --work-dir)               need_arg "$@"; WORK_ARG="$2"; shift 2 ;;
      --keep-build)             KEEP_BUILD=1; shift ;;
      --dry-run)                DRY_RUN=1; shift ;;
      --fortran-only)           FORTRAN_ONLY=1; shift ;;
      --netcdf-c)               need_arg "$@"; NETCDF_C="$2"; shift 2 ;;
      -h|--help)                usage 0 ;;
      *)                        printf 'unknown argument: %s\n' "$1" >&2; usage 1 ;;
   esac
done

[[ $DOWNLOAD_ONLY -eq 1 && $OFFLINE -eq 1 ]] && die "--download-only and --offline are mutually exclusive."
[[ $FORTRAN_ONLY -eq 1 && $USE_MPI -eq 1 ]] && die "--mpi applies to the netcdf-c build; it has no meaning with --fortran-only."
[[ "$JOBS" =~ ^[0-9]+$ ]] || die "--jobs must be a number, got '$JOBS'"

# Tarball checksums: explicit flag, else the built-in table.
[[ -n "$NC_SHA256" ]] || NC_SHA256="$(known_sha256 "netcdf-c-${NC_VERSION}")"
[[ -n "$NCF_SHA256" ]] || NCF_SHA256="$(known_sha256 "netcdf-fortran-${NCF_VERSION}")"

abspath() { (cd "$1" 2>/dev/null && pwd -P); }

# --- download / verify --------------------------------------------------
fetch() { # fetch URL DEST
   if command -v curl >/dev/null 2>&1; then
      curl -fsSL --retry 3 -o "$2" "$1"
   elif command -v wget >/dev/null 2>&1; then
      wget -q -O "$2" "$1"
   else
      die "neither curl nor wget is on PATH; fetch the tarballs elsewhere and use --offline."
   fi
}

sha256_of() {
   if command -v sha256sum >/dev/null 2>&1; then
      sha256sum "$1" | awk '{print $1}'
   else
      shasum -a 256 "$1" | awk '{print $1}'
   fi
}

# obtain NAME VERSION SHA DIR -> path of a verified tarball in DIR
obtain() {
   local name="$1" ver="$2" sha="$3" dir="$4"
   local file="${name}-${ver}.tar.gz"
   local url="${URL_BASE}/${name}/${ver}/${file}"
   local dest="${dir}/${file}"
   if [[ -f "$dest" ]]; then
      note "==> using ${dest}" >&2
   elif [[ $OFFLINE -eq 1 ]]; then
      die "--offline: ${dest} is missing.
Fetch it on a machine with internet access:
  tools/build_netcdf.sh --download-only --download-dir <dir> \\
      --netcdf-c-version ${NC_VERSION} --netcdf-fortran-version ${NCF_VERSION}"
   else
      note "==> downloading ${url}" >&2
      fetch "$url" "${dest}.part" || { rm -f "${dest}.part"; die "download failed: $url"; }
      mv "${dest}.part" "$dest"
   fi
   local got
   got="$(sha256_of "$dest")"
   if [[ -n "$sha" ]]; then
      [[ "$got" == "$sha" ]] || die "checksum mismatch for ${dest}
  expected  ${sha}
  got       ${got}
(delete the file to re-download it)"
      note "    sha256 OK" >&2
   else
      warn "${file} is not in the pinned checksum table; NOT verified.
         sha256 = ${got}   (pass it back with --${name}-sha256 to pin it)"
   fi
   printf '%s\n' "$dest"
}

if [[ $DOWNLOAD_ONLY -eq 1 ]]; then
   [[ -n "$DOWNLOAD_DIR" ]] || DOWNLOAD_DIR="$PWD/netcdf-src"
   mkdir -p "$DOWNLOAD_DIR"
   [[ $FORTRAN_ONLY -eq 1 ]] || obtain netcdf-c "$NC_VERSION" "$NC_SHA256" "$DOWNLOAD_DIR" >/dev/null
   obtain netcdf-fortran "$NCF_VERSION" "$NCF_SHA256" "$DOWNLOAD_DIR" >/dev/null
   note ""
   note "Tarballs are in $(abspath "$DOWNLOAD_DIR"). Build without network access with:"
   note "  tools/build_netcdf.sh --offline $(abspath "$DOWNLOAD_DIR") --fc <fc> --prefix <prefix> ..."
   exit 0
fi

# --- compilers ------------------------------------------------------------
# A Cray PE is recognised by its environment, not by the compiler name: the
# wrappers cc/ftn only exist with craype loaded, which sets CRAYPE_VERSION.
ON_CRAY=0
[[ -n "${CRAYPE_VERSION:-}" ]] && ON_CRAY=1

if [[ -z "$FC" && $ON_CRAY -eq 1 ]]; then FC=ftn; fi
if [[ $USE_MPI -eq 1 && $ON_CRAY -eq 0 && -z "$FC" ]]; then FC=mpifort; fi
[[ -n "$FC" ]] || die "no Fortran compiler given. Pass --fc <compiler> (or set FC).
Rakali supports gfortran (>= 15), nvfortran, ifx, and flang-new;
on a Cray PE use the ftn wrapper."
command -v "$FC" >/dev/null 2>&1 || die "Fortran compiler '$FC' is not on PATH.
Load your toolchain first (module load / source setvars.sh / conda activate)."
FC_PATH="$(command -v "$FC")"
FC_TAG="$(basename "$FC")"

# The C compiler paired with a Fortran compiler (netcdf-c is plain C, so any
# working C compiler is ABI-compatible; pairing just avoids surprises).
paired_cc() {
   local f
   f="$(basename "$1")"
   case "$f" in
      ftn)            echo cc ;;
      mpifort|mpif90) echo mpicc ;;
      gfortran-*)     echo "gcc-${f#gfortran-}" gcc ;;
      gfortran)       echo gcc ;;
      nvfortran|pgfortran|pgf90) echo nvc pgcc ;;
      ifx)            echo icx ;;
      ifort)          echo icc icx ;;
      flang-new*|flang*) echo clang ;;
      *)              echo cc gcc ;;
   esac
}

CC_PATH=""
if [[ -n "$CC_ARG" ]]; then
   command -v "$CC_ARG" >/dev/null 2>&1 || die "C compiler '$CC_ARG' is not on PATH."
   CC_PATH="$(command -v "$CC_ARG")"
elif [[ $FORTRAN_ONLY -eq 0 ]]; then
   if [[ -n "${CC:-}" ]] && command -v "$CC" >/dev/null 2>&1 \
         && { [[ $USE_MPI -eq 0 ]] || [[ $ON_CRAY -eq 1 ]]; }; then
      CC_PATH="$(command -v "$CC")"
   elif [[ $USE_MPI -eq 1 && $ON_CRAY -eq 0 ]]; then
      command -v mpicc >/dev/null 2>&1 || die "--mpi needs mpicc on PATH (or pass --cc)."
      CC_PATH="$(command -v mpicc)"
   else
      for c in $(paired_cc "$FC"); do
         if command -v "$c" >/dev/null 2>&1; then CC_PATH="$(command -v "$c")"; break; fi
      done
   fi
   [[ -n "$CC_PATH" ]] || die "no C compiler found for '$FC_TAG'. Pass --cc <compiler>."
fi
# --fortran-only never needed a C compiler of ours: CMake picks one unless
# --cc is given explicitly (keeps tools/build_netcdf_fortran.sh unchanged).

CC_TAG=""
[[ -n "$CC_PATH" ]] && CC_TAG="$(basename "$CC_PATH")"

# Capture --version output first: `cmd | grep -q` under pipefail can report
# a SIGPIPE'd cmd as a failure.
FC_BANNER="$("$FC_PATH" --version 2>/dev/null || true)"
IS_NVHPC=0
[[ "$FC_BANNER" =~ nvfortran|NVIDIA|nvidia|PGI ]] && IS_NVHPC=1

# Compiler flags. -fPIC: libtool does not know every compiler's PIC flag
# (nvc, cray clang behind cc), and a shared libnetcdf with non-PIC objects
# fails at link time. Every supported C compiler accepts -fPIC.
NC_CFLAGS="-O2 -fPIC"
# GCC >= 15 defaults to C23, in which `bool` is a keyword and `()` means
# "no arguments"; netcdf-c's code base predates that. Pin gnu17 (a
# precaution: 4.9.3 was only ever built here with it).
CC_BANNER=""
[[ -n "$CC_PATH" ]] && CC_BANNER="$("$CC_PATH" --version 2>/dev/null || true)"
if [[ "${CC_BANNER%%$'\n'*}" == *"(GCC)"* || "${CC_BANNER%%$'\n'*}" == gcc* ]]; then
   gcc_major="$("$CC_PATH" -dumpversion 2>/dev/null | cut -d. -f1)"
   if [[ "$gcc_major" =~ ^[0-9]+$ ]] && (( gcc_major >= 15 )); then
      NC_CFLAGS="$NC_CFLAGS -std=gnu17"
   fi
fi
NC_CFLAGS="$NC_CFLAGS${EXTRA_CFLAGS:+ $EXTRA_CFLAGS}"

# nvc/nvfortran target the BUILD host's CPU unless -tp says otherwise; a
# library built on a newer login node can SIGILL on older compute nodes.
if [[ $IS_NVHPC -eq 1 && "$EXTRA_CFLAGS $EXTRA_FFLAGS" != *-tp* ]]; then
   warn "NVHPC targets THIS host's CPU by default. If the compute nodes have a
         different CPU than this node, pass --cflags '-tp=<cpu>' --fflags '-tp=<cpu>'
         (and -DRAKALI_EXTRA_FORTRAN_FLAGS=-tp=<cpu> to Rakali)."
fi

# Cray static linking: older PEs (and CRAYPE_LINK_TYPE=static) link static.
if [[ $ON_CRAY -eq 1 && "${CRAYPE_LINK_TYPE:-dynamic}" == static && $STATIC -eq 0 ]]; then
   note "CRAYPE_LINK_TYPE=static -> building static libraries (--static)."
   STATIC=1
fi

# --- HDF5 (never built here) ----------------------------------------------
HDF5_INC=""
HDF5_LIB=""
HDF5_SRC=""

hdf5_from_prefix() { # prefix -> sets HDF5_INC/HDF5_LIB, returns 1 if not an HDF5
   local p="$1" d
   HDF5_INC=""; HDF5_LIB=""
   for d in "$p/include" "$p/include/hdf5/serial" "$p/include/hdf5/openmpi" \
            "$p/include/hdf5/mpich" "$p/include/hdf5"; do
      if [[ -f "$d/hdf5.h" ]]; then HDF5_INC="$d"; break; fi
   done
   for d in "$p/lib" "$p/lib64" "$p/lib/x86_64-linux-gnu/hdf5/serial" \
            "$p/lib/x86_64-linux-gnu" "$p/lib/aarch64-linux-gnu/hdf5/serial"; do
      if compgen -G "$d/libhdf5.so*" >/dev/null || [[ -f "$d/libhdf5.a" ]]; then
         HDF5_LIB="$d"; break
      fi
   done
   [[ -n "$HDF5_INC" && -n "$HDF5_LIB" ]]
}

hdf5_from_wrapper() { # h5cc|h5pcc -> parse -I/-L out of `-show`
   local w="$1" tok p
   command -v "$w" >/dev/null 2>&1 || return 1
   HDF5_INC=""; HDF5_LIB=""
   for tok in $("$w" -show 2>/dev/null); do
      case "$tok" in
         -I*) [[ -z "$HDF5_INC" && -f "${tok#-I}/hdf5.h" ]] && HDF5_INC="${tok#-I}" ;;
         -L*) [[ -z "$HDF5_LIB" ]] && HDF5_LIB="${tok#-L}" ;;
         */libhdf5.a|*/libhdf5.so) [[ -z "$HDF5_LIB" ]] && HDF5_LIB="$(dirname "$tok")" ;;
      esac
   done
   # Some wrappers omit -I for a default include dir; fall back to the
   # installation point `-showconfig` reports.
   if [[ -z "$HDF5_INC" || -z "$HDF5_LIB" ]]; then
      p="$("$w" -showconfig 2>/dev/null | awk -F': *' '/Installation point/{print $2; exit}')"
      if [[ -n "$p" ]]; then
         local inc="$HDF5_INC" lib="$HDF5_LIB"
         hdf5_from_prefix "$p" || true
         [[ -n "$inc" ]] && HDF5_INC="$inc"
         [[ -n "$lib" ]] && HDF5_LIB="$lib"
      fi
   fi
   [[ -n "$HDF5_INC" && -n "$HDF5_LIB" ]]
}

hdf5_from_pkgconfig() {
   local pc tok
   command -v pkg-config >/dev/null 2>&1 || return 1
   for pc in hdf5 hdf5-serial hdf5-openmpi hdf5-mpich; do
      pkg-config --exists "$pc" 2>/dev/null || continue
      HDF5_INC=""; HDF5_LIB=""
      for tok in $(pkg-config --cflags-only-I "$pc"); do
         [[ -z "$HDF5_INC" && -f "${tok#-I}/hdf5.h" ]] && HDF5_INC="${tok#-I}"
      done
      for tok in $(pkg-config --libs-only-L "$pc"); do
         [[ -z "$HDF5_LIB" ]] && HDF5_LIB="${tok#-L}"
      done
      if [[ -n "$HDF5_LIB" && -z "$HDF5_INC" && -f /usr/include/hdf5.h ]]; then HDF5_INC=/usr/include; fi
      if [[ -n "$HDF5_INC" && -z "$HDF5_LIB" ]]; then
         hdf5_from_prefix "$(pkg-config --variable=prefix "$pc")" || true
      fi
      [[ -n "$HDF5_INC" && -n "$HDF5_LIB" ]] && { HDF5_SRC="pkg-config ($pc)"; return 0; }
   done
   return 1
}

find_hdf5() {
   local w
   if [[ -n "$HDF5_ARG" ]]; then
      [[ -d "$HDF5_ARG" ]] || die "--hdf5 '$HDF5_ARG' does not exist."
      hdf5_from_prefix "$HDF5_ARG" || die "--hdf5 '$HDF5_ARG' has no hdf5.h + libhdf5 under include/ and lib/ (or lib64/)."
      HDF5_SRC="--hdf5"; return
   fi
   for v in HDF5_DIR HDF5_ROOT; do
      if [[ -n "${!v:-}" ]]; then
         hdf5_from_prefix "${!v}" || die "\$$v='${!v}' is set but has no hdf5.h + libhdf5 under it."
         HDF5_SRC="\$$v"; return
      fi
   done
   if [[ $USE_MPI -eq 1 ]]; then w="h5pcc"; else w="h5cc"; fi
   if hdf5_from_wrapper "$w"; then HDF5_SRC="$(command -v "$w") -show"; return; fi
   if hdf5_from_pkgconfig; then return; fi
   die "could not find HDF5, and this script does not build it.

Point it at an existing HDF5 install (serial is enough for Rakali):
  --hdf5 <prefix>                 prefix holding include/hdf5.h + lib/libhdf5*
  export HDF5_DIR=<prefix>        (set by Cray's cray-hdf5 module)
  module load cray-hdf5           NERSC / Cray PE
  module load hdf5                most HPC sites (name varies)
  sudo apt install libhdf5-dev    Debian/Ubuntu
  sudo dnf install hdf5-devel     Fedora/RHEL
Searched: --hdf5, \$HDF5_DIR, \$HDF5_ROOT, \`$w -show\`, pkg-config hdf5."
}

if [[ $FORTRAN_ONLY -eq 0 ]]; then
   find_hdf5
   HDF5_INC="$(abspath "$HDF5_INC")"
   HDF5_LIB="$(abspath "$HDF5_LIB")"
   [[ -f "$HDF5_INC/hdf5.h" ]] || die "HDF5 include dir '$HDF5_INC' has no hdf5.h."
   compgen -G "$HDF5_LIB/libhdf5_hl*" >/dev/null || die "HDF5 in '$HDF5_LIB' has no libhdf5_hl (the high-level library); netcdf-c requires it."
   HDF5_VERSION="$(awk '/define H5_VERSION/{gsub(/"/,"",$3); print $3; exit}' "$HDF5_INC/H5pubconf.h" 2>/dev/null || true)"
   HDF5_PARALLEL=0
   grep -qE '^#define H5_HAVE_PARALLEL 1' "$HDF5_INC/H5pubconf.h" 2>/dev/null && HDF5_PARALLEL=1
   if [[ $HDF5_PARALLEL -eq 1 && $USE_MPI -eq 0 ]]; then
      die "the HDF5 found via ${HDF5_SRC} (${HDF5_INC}) is PARALLEL.
A parallel HDF5 needs an MPI C compiler for netcdf-c. Either
  * pass --mpi (Rakali must then be built with -DRAKALI_ENABLE_MPI=ON), or
  * point at a serial HDF5 (e.g. module load cray-hdf5, not cray-hdf5-parallel)."
   fi
   if [[ $HDF5_PARALLEL -eq 0 && $USE_MPI -eq 1 ]]; then
      die "--mpi was given but the HDF5 found via ${HDF5_SRC} (${HDF5_INC}) is SERIAL.
Load a parallel HDF5 (cray-hdf5-parallel, hdf5-openmpi, ...) or drop --mpi --
Rakali only needs a serial netcdf."
   fi
fi

# --- netcdf-c for --fortran-only --------------------------------------------
if [[ $FORTRAN_ONLY -eq 1 ]]; then
   if [[ -z "$NETCDF_C" && -n "${NETCDFC_DIR:-}" ]]; then
      NETCDF_C="$NETCDFC_DIR"; NETCDF_C_SRC="\$NETCDFC_DIR"
   elif [[ -z "$NETCDF_C" ]] && command -v nc-config >/dev/null 2>&1; then
      NETCDF_C="$(nc-config --prefix)"; NETCDF_C_SRC="nc-config on PATH"
   else
      NETCDF_C_SRC="--netcdf-c"
   fi
   [[ -n "$NETCDF_C" ]] || die "could not find netcdf-c.

Give it explicitly with --netcdf-c <prefix>, or install one -- ANY C build
will do, it does not have to match your Fortran compiler:

  Debian/Ubuntu   sudo apt install libnetcdf-dev
  Fedora/RHEL     sudo dnf install netcdf-devel
  conda           conda install -c conda-forge libnetcdf
  HPC site        module load netcdf-c   (name varies)
  none at all     tools/build_netcdf.sh   (builds netcdf-c too)"
   [[ -d "$NETCDF_C" ]] || die "netcdf-c prefix '$NETCDF_C' does not exist (from $NETCDF_C_SRC)."
   [[ -f "$NETCDF_C/include/netcdf.h" ]] || die \
      "'$NETCDF_C' has no include/netcdf.h -- that is not a netcdf-c prefix (from $NETCDF_C_SRC)."
   # CMake runs from the build tree, where a relative prefix means nothing.
   NETCDF_C="$(abspath "$NETCDF_C")"
fi

# --- prefix ---------------------------------------------------------------
if [[ -z "$PREFIX" ]]; then
   if [[ $FORTRAN_ONLY -eq 1 ]]; then
      PREFIX="$HOME/opt/netcdf-fortran-${NCF_VERSION}-${FC_TAG}"
   else
      PREFIX="$HOME/opt/netcdf-${NC_VERSION}-${FC_TAG}"
   fi
fi
[[ "$PREFIX" == /* ]] || PREFIX="$PWD/$PREFIX"
if [[ $DRY_RUN -eq 0 ]]; then
   mkdir -p "$PREFIX"
   PREFIX="$(abspath "$PREFIX")"
fi
[[ $FORTRAN_ONLY -eq 1 ]] || NETCDF_C="$PREFIX"
NC_LIBDIR="$NETCDF_C/lib"
[[ -d "$NC_LIBDIR" || ! -d "$NETCDF_C/lib64" ]] || NC_LIBDIR="$NETCDF_C/lib64"

if [[ $STATIC -eq 1 ]]; then LINK_KIND="static"; else LINK_KIND="shared"; fi

note "-------------------------------------------------------------------"
[[ $FORTRAN_ONLY -eq 1 ]] || note "  netcdf-c         ${NC_VERSION}"
note "  netcdf-fortran   ${NCF_VERSION}"
note "  Fortran compiler ${FC_PATH}"
[[ -n "$CC_PATH" ]] && note "  C compiler       ${CC_PATH}"
[[ $ON_CRAY -eq 1 ]] && note "  Cray PE          yes (CRAYPE_VERSION=${CRAYPE_VERSION}, PE_ENV=${PE_ENV:-?})"
if [[ $FORTRAN_ONLY -eq 1 ]]; then
   note "  netcdf-c         ${NETCDF_C}  (found via ${NETCDF_C_SRC})"
else
   note "  HDF5             ${HDF5_VERSION:-unknown}  (found via ${HDF5_SRC})"
   note "                   include ${HDF5_INC}"
   note "                   lib     ${HDF5_LIB}"
   if [[ $USE_MPI -eq 1 ]]; then note "  parallel I/O     yes (--mpi)"; else note "  parallel I/O     no (serial)"; fi
fi
note "  libraries        ${LINK_KIND}"
note "  install prefix   ${PREFIX}"
note "  parallel jobs    ${JOBS}"
[[ $DRY_RUN -eq 1 ]] && note "  DRY RUN          nothing is downloaded, built or installed"
note "-------------------------------------------------------------------"

for tool in cmake make tar; do
   command -v "$tool" >/dev/null 2>&1 || die "'$tool' is required but not on PATH."
done

# --- work tree --------------------------------------------------------------
if [[ $DRY_RUN -eq 1 ]]; then
   WORK="${WORK_ARG:-<work-dir>}"
elif [[ -n "$WORK_ARG" ]]; then
   mkdir -p "$WORK_ARG"
   WORK="$(abspath "$WORK_ARG")"
else
   WORK="$(mktemp -d "${TMPDIR:-/tmp}/netcdf-build-XXXXXX")"
fi
cleanup() {
   if [[ $DRY_RUN -eq 1 ]]; then return; fi
   if [[ $KEEP_BUILD -eq 1 ]]; then
      note "build tree kept at ${WORK}"
   elif [[ -z "$WORK_ARG" ]]; then
      rm -rf "$WORK"
   else
      rm -rf "$WORK/netcdf-c-${NC_VERSION}" "$WORK/netcdf-fortran-${NCF_VERSION}" "$WORK/ncf-build" "$WORK/smoke"
   fi
}
trap cleanup EXIT

# indir DIR cmd... : run cmd in DIR without changing our own cwd.
indir() { (cd "$1" && shift && "$@"); }

# run LOGNAME cmd... : run quietly into $WORK/LOGNAME.log, show its tail on
# failure. Under --dry-run it only prints the command.
run() {
   local log="$1"; shift
   if [[ $DRY_RUN -eq 1 ]]; then
      printf '    $'; printf ' %q' "$@"; printf '\n'
      return 0
   fi
   "$@" > "$WORK/$log.log" 2>&1 || {
      tail -40 "$WORK/$log.log"
      die "$log failed; full log: $WORK/$log.log (re-run with --keep-build to keep it)"
   }
}

# --- sources ----------------------------------------------------------------
SRC_DIR="${DOWNLOAD_DIR:-$WORK}"
if [[ $DRY_RUN -eq 1 ]]; then
   [[ $FORTRAN_ONLY -eq 1 ]] || note "==> would fetch ${URL_BASE}/netcdf-c/${NC_VERSION}/netcdf-c-${NC_VERSION}.tar.gz (sha256 ${NC_SHA256:-UNPINNED})"
   note "==> would fetch ${URL_BASE}/netcdf-fortran/${NCF_VERSION}/netcdf-fortran-${NCF_VERSION}.tar.gz (sha256 ${NCF_SHA256:-UNPINNED})"
else
   mkdir -p "$SRC_DIR"
   if [[ $FORTRAN_ONLY -eq 0 ]]; then
      NC_TARBALL="$(obtain netcdf-c "$NC_VERSION" "$NC_SHA256" "$SRC_DIR")"
      tar xzf "$NC_TARBALL" -C "$WORK"
   fi
   NCF_TARBALL="$(obtain netcdf-fortran "$NCF_VERSION" "$NCF_SHA256" "$SRC_DIR")"
   tar xzf "$NCF_TARBALL" -C "$WORK"
fi

# --- netcdf-c (autotools) -----------------------------------------------------
if [[ $FORTRAN_ONLY -eq 0 ]]; then
   NC_SRC="$WORK/netcdf-c-${NC_VERSION}"
   [[ $DRY_RUN -eq 1 || -x "$NC_SRC/configure" ]] || die "unpacked netcdf-c tree has no configure at $NC_SRC"

   # Everything Rakali does not use is off: remote access (DAP, byte-range,
   # S3), NCZarr, filter plugins, libxml2, examples. That removes libcurl and
   # libxml2 from the dependency list. NetCDF-4/HDF5 stays on; szip is linked
   # only if the HDF5 itself was built with it (configure auto-detects).
   NC_CONF=(
      --prefix="$PREFIX"
      --enable-netcdf-4
      --disable-remote-functionality
      --disable-dap
      --disable-byterange
      --disable-nczarr
      --disable-libxml2
      --disable-plugins
      --disable-filter-testing
      --disable-filter-blosc
      --disable-filter-zstd
      --disable-examples
      --disable-dependency-tracking
   )
   [[ $CHECK -eq 1 ]] || NC_CONF+=(--disable-testsets)
   if [[ $STATIC -eq 1 ]]; then
      NC_CONF+=(--enable-static --disable-shared)
   else
      NC_CONF+=(--enable-shared --disable-static)
   fi
   [[ $USE_MPI -eq 1 ]] && NC_CONF+=(--enable-parallel4)

   # rpath: libnetcdf.so records where libhdf5 lives, so nothing downstream
   # (netcdf-fortran, Rakali, ncdump) needs LD_LIBRARY_PATH for it.
   NC_LDFLAGS="-L$HDF5_LIB"
   [[ $STATIC -eq 1 ]] || NC_LDFLAGS="$NC_LDFLAGS -Wl,-rpath,$HDF5_LIB -Wl,-rpath,$PREFIX/lib"

   note ""
   note "==> netcdf-c ${NC_VERSION}: configure"
   run nc-configure indir "$NC_SRC" env \
      CC="$CC_PATH" CFLAGS="$NC_CFLAGS" CPPFLAGS="-I$HDF5_INC" LDFLAGS="$NC_LDFLAGS" \
      ./configure "${NC_CONF[@]}"
   note "==> netcdf-c: build (-j${JOBS})"
   run nc-build make -C "$NC_SRC" -j "$JOBS"
   if [[ $CHECK -eq 1 ]]; then
      note "==> netcdf-c: make check"
      run nc-check make -C "$NC_SRC" -j "$JOBS" check
   fi
   note "==> netcdf-c: install"
   run nc-install make -C "$NC_SRC" install
   if [[ $DRY_RUN -eq 0 ]]; then
      [[ -f "$PREFIX/include/netcdf.h" ]] || die "netcdf-c install finished but $PREFIX/include/netcdf.h is missing."
      [[ "$("$PREFIX/bin/nc-config" --has-nc4)" == yes ]] || die "the installed netcdf-c has no NetCDF-4/HDF5 support."
   fi
fi

# --- netcdf-fortran (CMake) ---------------------------------------------------
NCF_SRC="$WORK/netcdf-fortran-${NCF_VERSION}"
NCF_BUILD="$WORK/ncf-build"
# netcdf-fortran finds netcdf-c through find_package(netCDF) first and falls
# back to FIND_LIBRARY; CMAKE_PREFIX_PATH satisfies both. The generated
# netCDF-FortranConfig.cmake then records THIS netcdf-c, which is what
# Rakali's find_package(netCDF-Fortran) resolves against.
NCF_CMAKE=(
   -S "$NCF_SRC" -B "$NCF_BUILD"
   -DCMAKE_Fortran_COMPILER="$FC_PATH"
   -DCMAKE_INSTALL_PREFIX="$PREFIX"
   -DCMAKE_INSTALL_LIBDIR=lib
   -DCMAKE_PREFIX_PATH="$NETCDF_C"
   -DCMAKE_BUILD_TYPE=Release
   -DBUILD_EXAMPLES=OFF
)
[[ -n "$CC_PATH" ]] && NCF_CMAKE+=(-DCMAKE_C_COMPILER="$CC_PATH")
[[ -n "$EXTRA_FFLAGS" ]] && NCF_CMAKE+=(-DCMAKE_Fortran_FLAGS="$EXTRA_FFLAGS")
[[ -n "$EXTRA_CFLAGS" ]] && NCF_CMAKE+=(-DCMAKE_C_FLAGS="$EXTRA_CFLAGS")
if [[ $CHECK -eq 1 ]]; then NCF_CMAKE+=(-DENABLE_TESTS=ON); else NCF_CMAKE+=(-DENABLE_TESTS=OFF); fi
if [[ $STATIC -eq 1 ]]; then
   NCF_CMAKE+=(-DBUILD_SHARED_LIBS=OFF)
else
   NCF_RPATH="$PREFIX/lib"
   [[ "$NC_LIBDIR" == "$PREFIX/lib" ]] || NCF_RPATH="$NCF_RPATH;$NC_LIBDIR"
   [[ -n "$HDF5_LIB" ]] && NCF_RPATH="$NCF_RPATH;$HDF5_LIB"
   NCF_CMAKE+=(-DBUILD_SHARED_LIBS=ON
      -DCMAKE_INSTALL_RPATH="$NCF_RPATH"
      -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON)
fi

note ""
note "==> netcdf-fortran ${NCF_VERSION}: configure"
run ncf-configure cmake "${NCF_CMAKE[@]}"
note "==> netcdf-fortran: build (-j${JOBS})"
run ncf-build cmake --build "$NCF_BUILD" -j "$JOBS"
if [[ $CHECK -eq 1 ]]; then
   note "==> netcdf-fortran: ctest"
   run ncf-check ctest --test-dir "$NCF_BUILD" --output-on-failure
fi
note "==> netcdf-fortran: install"
run ncf-install cmake --install "$NCF_BUILD"
[[ $DRY_RUN -eq 1 || -f "$PREFIX/include/netcdf.mod" ]] || die "install finished but $PREFIX/include/netcdf.mod is missing."

# --- prove it actually works before claiming success --------------------------
note "==> smoke test (compile, link, write a NetCDF-4 file, read it back)"
SMOKE_LIBS=(-L"$PREFIX/lib" -lnetcdff)
[[ "$NC_LIBDIR" == "$PREFIX/lib" ]] || SMOKE_LIBS+=(-L"$NC_LIBDIR")
SMOKE_LIBS+=(-lnetcdf)
# The executable gets an RPATH to the prefix, as CMake gives Rakali's;
# everything below it must then resolve through the libraries' own RPATHs.
if [[ $STATIC -eq 0 ]]; then
   SMOKE_LIBS+=(-Wl,-rpath,"$PREFIX/lib")
   [[ "$NC_LIBDIR" == "$PREFIX/lib" ]] || SMOKE_LIBS+=(-Wl,-rpath,"$NC_LIBDIR")
fi
if [[ $STATIC -eq 1 && -x "$NETCDF_C/bin/nc-config" && $DRY_RUN -eq 0 ]]; then
   # A static libnetcdf.a carries no dependency list; nc-config has it.
   read -r -a extra <<< "$("$NETCDF_C/bin/nc-config" --libs --static 2>/dev/null || "$NETCDF_C/bin/nc-config" --libs)"
   SMOKE_LIBS+=("${extra[@]}")
fi
if [[ $DRY_RUN -eq 0 ]]; then
   mkdir -p "$WORK/smoke"
   cat > "$WORK/smoke/smoke.f90" <<'EOF'
program smoke
   use netcdf
   implicit none
   integer :: ncid, dimid, varid, ierr
   real :: wrote(4) = [1.0, 2.0, 3.0, 4.0], readback(4)
   ierr = nf90_create("smoke.nc", ior(NF90_CLOBBER, NF90_NETCDF4), ncid); if (ierr /= nf90_noerr) stop 1
   ierr = nf90_def_dim(ncid, "x", 4, dimid);            if (ierr /= nf90_noerr) stop 2
   ierr = nf90_def_var(ncid, "v", NF90_FLOAT, [dimid], varid, deflate_level=1)
   if (ierr /= nf90_noerr) stop 6
   ierr = nf90_enddef(ncid)
   ierr = nf90_put_var(ncid, varid, wrote)
   ierr = nf90_close(ncid);                             if (ierr /= nf90_noerr) stop 3
   ierr = nf90_open("smoke.nc", NF90_NOWRITE, ncid);    if (ierr /= nf90_noerr) stop 4
   ierr = nf90_inq_varid(ncid, "v", varid)
   ierr = nf90_get_var(ncid, varid, readback)
   ierr = nf90_close(ncid)
   if (any(abs(readback - wrote) > 0.0)) stop 5
   print *, "netcdf-fortran smoke test OK, library ", trim(nf90_inq_libvers())
end program
EOF
fi
# The NetCDF-4 + deflate write above goes through HDF5, so a netcdf-c
# without working HDF5 support fails here rather than inside Rakali.
run smoke-compile indir "$WORK/smoke" \
   "$FC_PATH" -o smoke smoke.f90 -I"$PREFIX/include" "${SMOKE_LIBS[@]}"
RUN_LIBPATH="$PREFIX/lib"
[[ "$NC_LIBDIR" == "$PREFIX/lib" ]] || RUN_LIBPATH="$RUN_LIBPATH:$NC_LIBDIR"
[[ -z "$HDF5_LIB" ]] || RUN_LIBPATH="$RUN_LIBPATH:$HDF5_LIB"
if [[ $DRY_RUN -eq 0 ]]; then
   # No LD_LIBRARY_PATH here on purpose: the RPATHs must be enough.
   ( cd "$WORK/smoke" && ./smoke ) || {
      note "smoke test failed with RPATH only; retrying with LD_LIBRARY_PATH=$RUN_LIBPATH"
      ( cd "$WORK/smoke" && LD_LIBRARY_PATH="$RUN_LIBPATH:${LD_LIBRARY_PATH:-}" ./smoke ) \
         || die "smoke test compiled but failed at runtime."
      warn "the libraries only resolve with LD_LIBRARY_PATH set -- source the env file below."
   }
fi

# --- environment file + next steps -----------------------------------------
ENV_FILE="$PREFIX/rakali-netcdf.env"
if [[ $DRY_RUN -eq 0 ]]; then
   {
      printf '# Written by tools/build_netcdf.sh on %s\n' "$(date -u +%Y-%m-%dT%H:%MZ)"
      printf '# netcdf-c %s (%s), netcdf-fortran %s, Fortran compiler %s\n' \
         "$( [[ $FORTRAN_ONLY -eq 1 ]] && echo "at $NETCDF_C" || echo "$NC_VERSION")" \
         "$LINK_KIND" "$NCF_VERSION" "$FC_PATH"
      [[ -n "$HDF5_LIB" ]] && printf '# HDF5 %s from %s\n' "${HDF5_VERSION:-?}" "$HDF5_LIB"
      printf 'export NETCDF_DIR=%q\n' "$NETCDF_C"
      printf 'export NETCDFF_DIR=%q\n' "$PREFIX"
      printf 'export PATH=%q:"$PATH"\n' "$PREFIX/bin"
      printf 'export CMAKE_PREFIX_PATH=%q"${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"\n' \
         "$PREFIX$( [[ "$NETCDF_C" != "$PREFIX" ]] && echo ":$NETCDF_C")"
      printf 'export PKG_CONFIG_PATH=%q"${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"\n' \
         "$PREFIX/lib/pkgconfig$( [[ "$NETCDF_C" != "$PREFIX" ]] && echo ":$NC_LIBDIR/pkgconfig")"
      printf 'export LD_LIBRARY_PATH=%q"${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"\n' "$RUN_LIBPATH"
   } > "$ENV_FILE"
fi

CMAKE_PP="$PREFIX"
[[ "$NETCDF_C" != "$PREFIX" ]] && CMAKE_PP="$PREFIX;$NETCDF_C"

note ""
note "-------------------------------------------------------------------"
if [[ $DRY_RUN -eq 1 ]]; then
   note "dry run complete -- nothing was built."
else
   note "netcdf-fortran ${NCF_VERSION} installed for ${FC_TAG} in ${PREFIX}"
   [[ $FORTRAN_ONLY -eq 1 ]] || note "netcdf-c ${NC_VERSION} (${LINK_KIND}) installed in ${PREFIX}"
   note "environment file: ${ENV_FILE}"
fi
note ""
note "Configure Rakali against it (from the Rakali source tree):"
note ""
note "  source ${ENV_FILE}"
note "  cmake -B build_${FC_TAG} -S . \\"
note "        -DCMAKE_Fortran_COMPILER=${FC_PATH} \\"
[[ -n "$CC_PATH" ]] && note "        -DCMAKE_C_COMPILER=${CC_PATH} \\"
if [[ $USE_MPI -eq 1 ]]; then
   note "        -DRAKALI_ENABLE_MPI=ON \\"
fi
if [[ $STATIC -eq 1 ]]; then
   note "        -DCMAKE_Fortran_STANDARD_LIBRARIES=\"\$(${NETCDF_C}/bin/nc-config --libs --static)\" \\"
fi
note "        -DCMAKE_PREFIX_PATH='${CMAKE_PP}'"
note ""
note "Configure prints 'NetCDF-Fortran found via CMake config: <dir>' --"
note "check that <dir> is under ${PREFIX}."
if [[ $USE_MPI -eq 1 ]]; then
   note ""
   note "This netcdf-c links MPI (parallel HDF5), so Rakali must be an MPI"
   note "build (-DRAKALI_ENABLE_MPI=ON) with the same MPI."
fi
if [[ $STATIC -eq 0 ]]; then
   note ""
   if [[ $FORTRAN_ONLY -eq 1 ]]; then
      note "libnetcdff carries an RPATH to itself and to ${NC_LIBDIR}; whether"
      note "that netcdf-c finds its HDF5 without LD_LIBRARY_PATH is up to how it"
      note "was built. The env file sets LD_LIBRARY_PATH if you need it."
   else
      note "The libraries carry an RPATH to themselves and to HDF5, so a built"
      note "executable runs without LD_LIBRARY_PATH. If you relocate a prefix,"
      note "source the env file (it sets LD_LIBRARY_PATH)."
   fi
fi
note "-------------------------------------------------------------------"
