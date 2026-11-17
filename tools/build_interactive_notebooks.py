#!/usr/bin/env python
"""Builds the interactive notebook(s) from tools/lib (code) + tools/interactive_text.py (narrative).

    python tools/build_interactive_notebooks.py          → CytokineDictionary_Practical_v3_interactive.ipynb (English)

The .ipynb files are GENERATED: edit the sources in tools/, never the notebooks."""
import sys
from pathlib import Path

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from interactive_cells import EXERCISES, PARTS, TITLES_ES, exercise_cell_source, setup_cell_source
from interactive_text_es import ES
from interactive_text import BADGE, COLAB_URL, TEXT

TEXT = {**TEXT, "es": ES}
NOTEBOOKS = {"en": "CytokineDictionary_Practical_v3_interactive.ipynb", "es": "CytokineDictionary_Practical_v3_interactive_es.ipynb"}
HIDDEN = {"cellView": "form", "jupyter": {"source_hidden": True}}


def build(lang="en"):
    T, name = TEXT[lang], NOTEBOOKS[lang]
    badge = BADGE.format(url=COLAB_URL.format(notebook=name))
    cells = [new_markdown_cell(T["title"].format(badge=badge)), new_markdown_cell(T["setup_md"]),
             new_code_cell(setup_cell_source(lang), metadata=dict(HIDDEN))]
    for i, (key, title, function, try_key) in enumerate(EXERCISES):
        title = TITLES_ES[i] if lang == "es" else title
        if key in PARTS:
            cells.append(new_markdown_cell(T[PARTS[key]]))
        cells.append(new_markdown_cell(T[key]))
        cells.append(new_code_cell(exercise_cell_source(title, function, lang), metadata=dict(HIDDEN)))
        cells.append(new_markdown_cell(T[try_key]))
    cells.append(new_markdown_cell(T["wrapup"]))
    nb = new_notebook(cells=cells, metadata={
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
        "colab": {"provenance": [], "toc_visible": True, "name": name},
    })
    nbformat.validate(nb)
    out = ROOT / name
    nbformat.write(nb, out)
    print(f"✅ {out.name}: {len(cells)} cells, {out.stat().st_size / 1e3:.0f} kB")
    return out


if __name__ == "__main__":
    for lang in (sys.argv[1:] or ["en", "es"]):
        build(lang)
