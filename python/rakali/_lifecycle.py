"""Shared lifecycle mixin for the Rakali simulation wrappers.

Both the structured (:class:`rakali.Rakali`) and unstructured
(:class:`rakali.unstr.UnstructuredSimulation`) front ends are opaque-handle
FFI wrappers with the same teardown shape: release the C handle exactly
once, support ``with`` blocks, and free on garbage collection. That common
shape lives here so the two classes feel identical and don't drift.

A subclass implements :meth:`_destroy_handle` to free its backend handle.
That method MUST be idempotent (guard on the subclass's own ``_alive`` /
``_destroyed`` flag) and safe to call on a half-built object — every public
teardown path funnels through it.
"""


def _owned_view(address, ctype, count, owner):
    """Zero-copy numpy view of ``count`` ``ctype`` values at ``address``.

    The view holds a reference to ``owner`` (the simulation wrapper), so the
    wrapper's ``__del__`` -- which frees the Fortran state -- cannot run while
    any view is alive. Calling ``close()`` explicitly still invalidates all
    views, as documented.
    """
    import numpy as np

    buf = (ctype * count).from_address(address)
    buf._owner = owner  # kept alive by numpy through the buffer base chain
    return np.frombuffer(buf, dtype=np.dtype(ctype))


def _require_file(path):
    """Raise FileNotFoundError if ``path`` is not an existing file.

    NetCDF open failures end in a Fortran ERROR STOP that kills the
    interpreter, so every path handed to the library is checked first.
    """
    import os

    if not os.path.isfile(str(path)):
        raise FileNotFoundError(f"File not found: {path}")


def _require_writable_dir(path):
    """Raise if the directory that would hold ``path`` is missing/unwritable."""
    import os

    parent = os.path.dirname(os.path.abspath(str(path)))
    if not os.path.isdir(parent):
        raise FileNotFoundError(f"Directory does not exist: {parent}")
    if not os.access(parent, os.W_OK):
        raise PermissionError(f"Directory is not writable: {parent}")


class _SimBase:
    """Context-manager + deterministic teardown for handle wrappers."""

    def _destroy_handle(self):
        """Free the backend C handle (idempotent). Subclass-provided."""
        raise NotImplementedError

    def close(self):
        """Release the solver handle. Idempotent; safe to call repeatedly."""
        self._destroy_handle()

    def destroy(self):
        """Alias for :meth:`close` (backwards-compatible name)."""
        self._destroy_handle()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self._destroy_handle()
        return False

    def __del__(self):
        # _destroy_handle guards on the subclass flag, so a partially built
        # object (flag unset) tears down cleanly; swallow shutdown races.
        try:
            self._destroy_handle()
        except Exception:
            pass
