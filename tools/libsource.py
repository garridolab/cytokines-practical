"""Assembles the helper library (the files of tools/lib, in alphabetical order) into ONE source text.

The same text is (1) pasted into the hidden setup cell of the notebook, (2) executed by the pre-computation script and
(3) executed by the smoke test, so the notebook, the cache and the tests always run exactly the same code."""
from pathlib import Path

LIB_DIR = Path(__file__).resolve().parent / "lib"


def library_source(exclude=()) -> str:
    """exclude: file-name prefixes to leave out (e.g. ("40_", "41_") for the notebook without the differential-expression exercises)."""
    parts = []
    for p in sorted(LIB_DIR.glob("*.py")):
        if p.name.startswith(tuple(exclude)):
            continue
        parts.append(p.read_text(encoding="utf-8").rstrip() + "\n")
    return "\n\n".join(parts)


def load_library(namespace=None) -> dict:
    """Execute the library in `namespace` (a fresh dict by default) and return it."""
    ns = {"__name__": "practical_lib"} if namespace is None else namespace
    exec(compile(library_source(), "practical_lib.py", "exec"), ns)
    return ns
