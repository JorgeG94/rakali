#!/usr/bin/env python3
"""Drift gate for the memory-accounting `bytes()` functions.

Fortran has no reflection, so a per-slot `bytes()` count is only as complete
as the human who wrote it: add an allocatable array (or a nested slot) to a
type and forget the matching term, and the counted footprint silently
under-reports.  The runtime reconciliation in `rki_mem_report`
(counted-vs-measured device mapping) catches it eventually — but only on a
GPU run, and only once.  This hook catches it at commit time.

RULE 1 — completeness (opt-in per type; a type is only checked once it HAS a
`bytes()`): for every derived type that binds `procedure :: bytes => <fn>`,
every COUNTABLE component must appear as `this%<name>` somewhere in <fn>.  A
component is countable when it is

  * an intrinsic allocatable (`real/integer/logical/complex, allocatable`), or
  * a component whose own declared type itself has a `bytes()` (nested slots
    and shared buffers — e.g. `scratch_3d_buffer_t`, `tracer_t`).

Derived components whose type has NO `bytes()` are not required by rule 1 —
rule 2 is what stops that being a free pass.

RULE 2 — existence (the blind spot rule 1 leaves): a type that has NO
`bytes()` at all used to be INVISIBLE to this hook, however much device
memory it mapped.  That is how `diag_var_t` (~3.2 GB), `bt_wide_t`
(~436 MB), `evp_workspace_t` (~154 MB) and `diag_mask_t` all went uncounted
while the hook stayed green.  So: if a type declares an intrinsic
allocatable component AND that component name appears in an `!$acc enter
data` / `!$acc declare` clause in the same file, the type is DEVICE-MAPPED
and must have a `bytes()` — or be listed in `DEVICE_MAPPED_EXEMPT` below
with a reason.

`DEVICE_MAPPED_EXEMPT` is a frozen BASELINE of the types that were already
device-mapped-and-uncounted when rule 2 was introduced.  It is not a
blessing — it is the to-do list.  Do not add to it to silence a new type;
give the new type a `bytes()` and wire the term into its parent's total.

Still NOT checked (known remaining gap): module-level `save` allocatables.
They have no owning type, so neither rule reaches them — every device
buffer in `src/comm/**` and `tracer_scratch_t`'s module storage is outside
this hook.  Only the runtime counted-vs-measured reconciliation in
`rki_mem_report` sees those.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"

DEVICE_MAPPED_EXEMPT = {
    # Pre-existing at the time rule 2 landed.  Each is device-mapped and
    # carries no `bytes()`; the reason records what is known, not that the
    # gap is acceptable.
    "data_input_field_t":
        "counted by hand — the reader walks fields(i)%f0/f1 itself",
    "vcoord_t":
        "coastal ALE vcoord slot; baseline gap (see docs/memaudit) — not on the "
        "state bytes() walk",
    "mesh_unstr_t":
        "unstructured backend; baseline gap — the unstructured state has no "
        "bytes() walk at all",
    "unstr_barotropic_state_t":
        "unstructured backend; baseline gap",
    "unstr_surface_forcing_t":
        "unstructured backend; baseline gap",
    "unstr_multilayer_state_t":
        "unstructured backend; baseline gap",
    "semi_implicit_state_unstr_t":
        "unstructured semi-implicit; baseline gap",
    "mg_level_t":
        "structured semi-implicit multigrid hierarchy; baseline gap "
        "(already named as a documented under-count by rule 1)",
    "agglom_level_unstr_t":
        "unstructured semi-implicit multigrid hierarchy; baseline gap",
}

RE_ACC_DIRECTIVE = re.compile(r"^\s*!\$acc\s+(?:enter\s+data|declare)\b", re.IGNORECASE)

RE_TYPE_DEF = re.compile(r"^\s*type\b(?!\s*\()(?:\s*,\s*[^:]+?)?\s*::\s*(\w+)", re.IGNORECASE)
RE_END_TYPE = re.compile(r"^\s*end\s*type\b", re.IGNORECASE)
RE_INTRINSIC = re.compile(
    r"^\s*(?:real|integer|logical|complex)\b[^:]*,\s*allocatable\s*::\s*(.+)$", re.IGNORECASE)
RE_DERIVED = re.compile(r"^\s*type\s*\(\s*(\w+)\s*\)[^:]*::\s*(.+)$", re.IGNORECASE)
RE_BYTES_BIND = re.compile(r"^\s*procedure\b.*::\s*bytes\s*=>\s*(\w+)", re.IGNORECASE)
RE_FUNC = re.compile(
    r"function\s+(\w+)\s*\(.*?\bend\s+function\s+\1\b", re.IGNORECASE | re.DOTALL)


def names_of(rhs):
    """Component names on a declaration RHS (strip dims + trailing comment)."""
    rhs = rhs.split("!", 1)[0]
    rhs = re.sub(r"\([^()]*\)", "", rhs)          # drop (:,:,:) / (0:nz) dims
    return [n.strip() for n in rhs.split(",") if n.strip()]


def acc_mapped_components(lines):
    """Component names appearing in any `!$acc enter data` / `declare` clause.

    Continuation lines are folded in, so a multi-line `enter data copyin(a, &
    b)` contributes both names.  Returned lower-cased.
    """
    mapped, i, n = set(), 0, len(lines)
    while i < n:
        if RE_ACC_DIRECTIVE.match(lines[i]):
            buf = lines[i]
            while buf.rstrip().endswith("&") and i + 1 < n:
                i += 1
                buf += " " + lines[i]
            mapped.update(m.group(1).lower() for m in re.finditer(r"%\s*(\w+)", buf))
        i += 1
    return mapped


def parse():
    """Return (types_with_bytes{type->fn}, components{type->[(name,dtype|None)]},
    func_bodies{fn->text}, device_mapped{type->(component, path)})."""
    types_with_bytes, components, func_bodies = {}, {}, {}
    device_mapped = {}
    for path in SRC.rglob("*.F90"):
        text = path.read_text()
        for m in RE_FUNC.finditer(text):
            func_bodies[m.group(1).lower()] = m.group(0)
        lines = text.splitlines()
        acc_names = acc_mapped_components(lines)
        cur = None
        for line in lines:
            if cur is None:
                md = RE_TYPE_DEF.match(line)
                if md:
                    cur = md.group(1)
                    components.setdefault(cur, [])
                continue
            if RE_END_TYPE.match(line):
                cur = None
                continue
            mb = RE_BYTES_BIND.match(line)
            if mb:
                types_with_bytes[cur] = mb.group(1).lower()
                continue
            mi = RE_INTRINSIC.match(line)
            if mi:
                for n in names_of(mi.group(1)):
                    components[cur].append((n, None))
                    if n.lower() in acc_names and cur not in device_mapped:
                        device_mapped[cur] = (n, str(path.relative_to(ROOT)))
                continue
            mdv = RE_DERIVED.match(line)
            if mdv:
                dtype = mdv.group(1)
                for n in names_of(mdv.group(2)):
                    components[cur].append((n, dtype))
    return types_with_bytes, components, func_bodies, device_mapped


def main():
    types_with_bytes, components, func_bodies, device_mapped = parse()
    have_bytes = {t.lower() for t in types_with_bytes}
    errors = []
    missing_bytes = []

    # ---- RULE 2: a device-mapped type must HAVE a bytes() ----------------
    for typ, (comp, path) in sorted(device_mapped.items()):
        if typ.lower() in have_bytes:
            continue
        if typ.lower() in {k.lower() for k in DEVICE_MAPPED_EXEMPT}:
            continue
        missing_bytes.append(
            f"{typ} ({path}): allocatable component `{comp}` is `!$acc enter "
            f"data`-mapped but the type has NO bytes() — device memory no "
            f"count can see. Add `procedure :: bytes => <fn>` and wire the "
            f"term into the parent slot's total.")

    # ---- RULE 1: a bytes() must be COMPLETE ------------------------------
    for typ, fn in sorted(types_with_bytes.items()):
        body = func_bodies.get(fn)
        if body is None:
            errors.append(f"{typ}: bytes() binds `{fn}` but no such function was found")
            continue
        for name, dtype in components.get(typ, []):
            countable = dtype is None or (dtype.lower() in have_bytes)
            if not countable:
                continue
            if not re.search(r"this\s*%\s*" + re.escape(name) + r"\b", body, re.IGNORECASE):
                kind = "array" if dtype is None else f"nested {dtype} slot"
                errors.append(
                    f"{typ}: {kind} `{name}` is not counted in {fn}() "
                    f"— add `+ arr_bytes(this%{name})`"
                    + ("" if dtype is None else f" or `+ this%{name}%bytes()`"))

    rc = 0
    if missing_bytes:
        print("Memory-accounting blind spot — a device-mapped type has no bytes():\n")
        for e in missing_bytes:
            print("  " + e)
        print("\n(DEVICE_MAPPED_EXEMPT in this file is a frozen baseline of "
              "pre-existing gaps, not a place to add new ones.)")
        rc = 1
    if errors:
        if missing_bytes:
            print()
        print("Memory-accounting drift — a countable component has no bytes() term:\n")
        for e in errors:
            print("  " + e)
        print("\n(If a component is intentionally uncounted, it must be a derived type "
              "without its own bytes(); intrinsic arrays are always required.)")
        rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
