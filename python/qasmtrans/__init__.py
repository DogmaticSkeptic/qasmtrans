"""Python package shim for qasmtrans."""
from importlib import import_module
import importlib.resources as resources
import sys
from pathlib import Path


def _import_core():
    # Try standard import first
    try:
        return import_module("qasmtrans_core")
    except ModuleNotFoundError:
        # Fallback: look for bundled binary in package resources
        try:
            with resources.path(__package__ or "qasmtrans", "qasmtrans_core") as p:
                sys.path.insert(0, str(Path(p).parent))
                return import_module("qasmtrans_core")
        except Exception:
            raise


_core = _import_core()

transpile_qasm = _core.transpile_qasm
emit_qick = _core.emit_qick
load_qick_config = _core.load_qick_config
TranspileOptions = _core.TranspileOptions
TranspileResult = _core.TranspileResult

__all__ = [
    "transpile_qasm",
    "emit_qick",
    "load_qick_config",
    "TranspileOptions",
    "TranspileResult",
]
