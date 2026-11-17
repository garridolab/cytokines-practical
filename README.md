# cytokines-practical

A ~2.5 h university practical on **transcriptomics data analysis** with the single-cell cytokine dictionary
(Oesinghaus *et al.* 2025: human PBMCs stimulated with 90 cytokines), pseudobulked per donor × cytokine × cell type.

| Version | Notebook | Slides |
|---|---|---|
| Coding practical (EN / ES) | `CytokineDictionary_Practical_v3.ipynb`, `…_v3_es.ipynb` | `slides/v3`, `slides/v3_es` |
| **Interactive, no-code practical (EN)** | **`CytokineDictionary_Practical_v3_interactive.ipynb`** | **`slides/v3_interactive`** (introduction + quiz QR only) |
| **Interactive, no-code practical (ES)** | **`CytokineDictionary_Practical_v3_interactive_es.ipynb`** | **`slides/v3_interactive_es`** |

## The interactive notebook

Students never write code: every exercise is **one grey cell** (code hidden) that shows a small panel of **controls**
(menus, sliders, text boxes, lists) and **figures/tables that redraw by themselves** whenever a control changes.
Each exercise has a short explanation, sensible defaults that already show a meaningful result, and 2–4 questions with folded answers.

| Part | Exercises |
|---|---|
| A · Python warm-up | 1 Variables and `print` · 2 Basic plots (simulated data) |
| B · Meet the data | 3 The AnnData object (`X`, `obs`, `var`, `layers`; build your own subset) · 4 Distributions of raw counts (boxplots) |
| C · Normalization | 5 Raw counts vs CPM · 6 The logarithm (base, pseudocount, variance) |
| D · Patterns | 7 Heatmaps (raw vs normalized) · 8 PCA: rotation and projection with two genes · 9 PCA of all samples · 10 PCA within one cell type |
| E · Differential expression | 11 Volcano plots · 12 Comparing two responses (correlation of log₂ fold changes) |

Differential expression is **never computed by the students**: the 87 cytokines × 4 cell types (348 DESeq2 comparisons, design `~ donor + cytokine`) were
computed once, **on ~500 immune genes only** (the 506 expressed genes of the cytokine-signalling, T-cell-activation and interferon gene sets, as in v3), and stored in
`cache/de_results.npz` (0.6 MB; it takes ~15 s to recompute).

### Run it

**Google Colab (what the students use).** Open the notebook from GitHub (badge in its first cell), press ▶ on the *Setup* cell and wait for `✅ Ready`
(~1–2 min: it installs the missing packages, downloads the 500 MB `.h5ad` from Google Drive and the small `de_results.npz` from this repository, normalizes and runs the PCAs).
Then press ▶ once in every exercise.

> ⚠️ **`cache/de_results.npz` must be committed and pushed to `main`** (it is whitelisted in `.gitignore`): Colab downloads it from
> `https://raw.githubusercontent.com/garridolab/cytokines-practical/main/cache/de_results.npz`. If the file cannot be downloaded the notebook still works:
> every comparison is then computed on demand with pydeseq2 (≈ 1 s each, pydeseq2 is installed automatically).

**Locally** (JupyterLab ≥ 3, Notebook ≥ 7 or VS Code):

```bash
pip install jupyterlab ipywidgets numpy pandas scipy matplotlib anndata gdown adjustText     # pydeseq2 is only needed without the cache
jupyter lab CytokineDictionary_Practical_v3_interactive.ipynb
```

Put `pseudobulks_4ct_proteincoding.h5ad` and `cache/de_results.npz` next to the notebook (or let the Setup cell download them). Tested with ipywidgets 7.7.1 (the version of Colab) and 8.1,
matplotlib 3.9–3.11, pandas 2.2–3.1, numpy 2.0–2.5, anndata 0.12–0.13.

### Develop it

The `.ipynb` is **generated**: do not edit it by hand.

```
tools/lib/*.py                      the Python helpers, one file per exercise (concatenated, in alphabetical order, into the hidden Setup cell)
tools/interactive_text.py           the narrative: explanations, questions and answers (EN); interactive_text_es.py = Spanish
tools/translations_es.py            Spanish version of every string of tools/lib (applied to the library SOURCE at build time by translate_lib.py)
tools/interactive_cells.py          the cells (install snippet, one-line exercise cells, order of the exercises)
tools/build_interactive_notebooks.py   → CytokineDictionary_Practical_v3_interactive.ipynb
tools/precompute_de.py              → cache/de_results.npz   (DESeq2 on the immune genes for every cytokine × cell type; ~15 s)
tools/smoke_test_interactive.py     re-calls every live render function with variations of every control; prints `ERRORS: 0`
```

```bash
python tools/build_interactive_notebooks.py [en|es]              # rebuild the notebook(s)
python tools/smoke_test_interactive.py figs --execute [--es]     # test (+ run the notebook from a fresh kernel); look at the PNGs in figs/
python tools/smoke_test_interactive.py figs --fallback           # also check the on-demand DESeq2 against the cache
python tools/precompute_de.py --workers 16                       # only if the data or the DE settings change
cd slides/v3_interactive[_es] && latexmk -xelatex main.tex       # slides
```

The numbers quoted in the answers of the questions (e.g. "IFN-γ in monocytes: 96 up / 95 down", "IL-4 vs IL-13: r = 0.97") were checked against the data;
the smoke test re-checks the main ones. Re-check them if you change a default or the data.
