#!/usr/bin/env python
"""Smoke test of the interactive notebook.

    python tools/smoke_test_interactive.py [figs_dir] [--execute] [--fallback] [--es]

1. Loads the helper library (the very same code that goes into the notebook's setup cell) and the real data.
2. Builds every exercise and re-calls its live render function with variations of EVERY control (all the options of the small
   menus, a sample of the long ones, minimum / middle / maximum of the sliders, empty / wrong / odd text, empty selections…).
3. Saves a PNG of every figure into `figs_dir` (default: tools/_smoke_figs) so that they can be looked at.
4. `--execute`: also runs the generated .ipynb from a fresh kernel with nbclient (checks that all cells run, and that no exercise
   recorded an error while drawing its initial state).
5. `--fallback`: also checks that the on-demand DESeq2 (used if the cache file is missing) reproduces the cached results.
Prints `ERRORS: 0` when everything is fine. Needs: ipywidgets, pydeseq2 (only for --fallback), nbclient (only for --execute).
"""
import itertools, os, sys, time, traceback
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from libsource import load_library

LANG = "es" if "--es" in sys.argv else "en"
args = [a for a in sys.argv[1:] if not a.startswith("--")]
FIGDIR = Path(args[0]) if args else HERE / "_smoke_figs"
FIGDIR.mkdir(parents=True, exist_ok=True)

if LANG == "es":
    from interactive_cells import lib_source
    L = {"__name__": "practical_lib"}
    exec(compile(lib_source("es"), "practical_lib_es.py", "exec"), L)
else:
    L = load_library()
import ipywidgets
import matplotlib
print(f"ipywidgets {ipywidgets.__version__} · matplotlib {matplotlib.__version__} · pandas {L['pd'].__version__} · numpy {L['np'].__version__} · anndata {L['ad'].__version__}")
L["display"] = lambda *a, **k: None                       # nothing to display without a notebook frontend
L["STRICT"] = True                                        # raise the errors instead of showing a friendly message
state = {"name": "x", "i": 0, "figs": 0}


def sink(fig):
    state["figs"] += 1
    if state["i"] < 3:                                    # keep the first 3 figures of each exercise
        fig.savefig(FIGDIR / f"{state['name']}_{state['i']}.png", dpi=70, bbox_inches="tight")
        state["i"] += 1


L["FIG_SINK"] = sink
t0 = time.time()
L["init_practical"](verbose=False, with_de=True)
adata, np = L["adata"], L["np"]
print(f"setup: {time.time() - t0:.1f} s · {adata.shape[0]:,} samples × {adata.shape[1]:,} genes · cached DE: {L['DE'].cached if L['DE'] is not None else 'not used (ES)'}")

errors, n_calls = [], 0
rng = np.random.default_rng(0)
BAD_TEXT = ["", "   ", "nonsense", "stat1", "STAT1,STAT1, ifit3 ; CXCL9\nIDO1", "ÄÖÜ ñ 日本語", "<b>x</b>", "A" * 300]


def variations(ctrl):
    """Values to try for one control (besides its default)."""
    W = ipywidgets
    if isinstance(ctrl, (W.Dropdown, W.ToggleButtons, W.SelectionSlider)):
        vals = list(ctrl._options_values)
        return vals if len(vals) <= 9 else [vals[0], vals[-1]] + [vals[i] for i in rng.choice(len(vals), 6, replace=False)]
    if isinstance(ctrl, W.SelectMultiple):
        vals = list(ctrl._options_values)
        return [(), (vals[0],), tuple(vals[:3]), tuple(vals[::max(1, len(vals) // 12)]), tuple(vals)]
    if isinstance(ctrl, (W.IntSlider, W.FloatSlider)):
        return [ctrl.min, (ctrl.min + ctrl.max) / 2 if isinstance(ctrl, W.FloatSlider) else int((ctrl.min + ctrl.max) / 2), ctrl.max]
    if isinstance(ctrl, W.Checkbox):
        return [True, False]
    if isinstance(ctrl, (W.Text, W.Combobox, W.Textarea)):
        return BAD_TEXT
    if isinstance(ctrl, W.IntText):
        return [0, 7, 123]
    return []


def call(ex, name, values):
    global n_calls
    n_calls += 1
    state["name"], state["i"] = name, 0
    try:
        ex.render(**values)
    except Exception:
        errors.append((name, values, traceback.format_exc(limit=4)))


def test_exercise(name, factory, extra=None):
    ex = factory()
    base = ex.values()
    t = time.time()
    call(ex, f"{name}_default", dict(base))
    for key, ctrl in ex.controls.items():
        if key.startswith("_"):
            continue
        for k, v in enumerate(variations(ctrl)):
            vals = dict(base)
            if key == "samples" and isinstance(v, str):
                pass
            vals[key] = v
            call(ex, f"{name}_{key}_{k}", vals)
    if extra:
        extra(ex)
    print(f"  {name:28s} {time.time() - t:5.1f} s")
    return ex


print("\nExercises:")
test_exercise("ex01_variables", L["exercise_variables"])
test_exercise("ex02_plotting", L["exercise_plotting"])


def ex3_extra(ex):
    sp = ex.picker
    base = ex.values()
    sp.f_cyt.value = "IL-4"
    sp.f_donor.value = "Donor3"
    assert len(sp.select.options) <= sp.max_options + len(sp.chosen)
    assert set(sp.chosen) <= set(sp.select.options), "the selection must survive the filters"
    sp.select.value = tuple(sp.select.options[:3])
    assert ex.controls["samples"].value.count("\n") == 2, "the hidden selection must follow the list"
    for n in (1, 2, 40, 200):
        ids = list(rng.choice(adata.obs_names.to_numpy(), n, replace=False))
        call(ex, f"ex03_n{n}", {**base, "samples": "\n".join(ids)})
    call(ex, "ex03_badid", {**base, "samples": "Donor99|X|Y\nDonor1|CD4_T_cell|PBS"})
    sp.f_cyt.value = sp.f_donor.value = sp.ALL


test_exercise("ex03_anndata", L["exercise_anndata"], ex3_extra)
test_exercise("ex04_boxplots", L["exercise_boxplots"])
test_exercise("ex05_normalization", L["exercise_normalization"])
test_exercise("ex06_logtransform", L["exercise_logtransform"])
test_exercise("ex07_heatmaps", L["exercise_heatmaps"])
test_exercise("ex08_pca_geometry", L["exercise_pca_geometry"])
test_exercise("ex09_pca_all", L["exercise_pca_all"])
test_exercise("ex10_pca_celltype", L["exercise_pca_celltype"])
if True:
    test_exercise("ex11_volcano", L["exercise_volcano"])


def ex12_extra(ex):
    for label, v in L["EXAMPLES12"]:
        if v:
            ex.set_values(cell_a=v[0], cyt_a=v[1], cell_b=v[2], cyt_b=v[3])
            call(ex, "ex12_example", ex.values())
    ex.set_values(cell_a="CD14_Mono", cyt_a="IL-4", cell_b="CD14_Mono", cyt_b="IL-13")
    # the correlation of the headline examples must match the numbers quoted in the narrative
    D = L["DE"]
    for (a, b), expected in {(("CD14_Mono", "IL-4"), ("CD14_Mono", "IL-13")): 0.97, (("CD4_T_cell", "IL-4"), ("CD4_T_cell", "IL-13")): 0.41,
                             (("CD14_Mono", "IFN-beta"), ("CD14_Mono", "IFN-gamma")): 0.78, (("CD14_Mono", "IFN-gamma"), ("CD14_Mono", "IL-4")): -0.35}.items():
        da, db = L["de_summary"](D.get(*a), 0.05, 0.5), L["de_summary"](D.get(*b), 0.05, 0.5)
        m, use = L["compare_tables"](da, db, "all", 0.05, 0.5)
        r = L["correlations"](m, use)[0]
        if abs(r - expected) > 0.006:
            errors.append(("narrative numbers", (a, b), f"r = {r:.3f}, the text says {expected}"))


if True:
    test_exercise("ex12_compare", L["exercise_compare"], ex12_extra)

# ── numbers quoted in the narrative (re-check when data or defaults change) ───────────────────────────────
print("\nNumbers quoted in the text:")
D = L["DE"]
checks = {}
if D is not None:
    d = L["de_summary"](D.get("CD14_Mono", "IFN-gamma"), 0.05, 0.5)
    checks = {"IFN-γ monocytes up": ((d["direction"] == "up").sum(), 96), "IFN-γ monocytes down": ((d["direction"] == "down").sum(), 95),
              "IFN-γ monocytes tested": (d["tested"].sum(), 452)}
    for thr, exp in ((0, 311), (1, 70), (2, 14)):
        dd = L["de_summary"](D.get("CD14_Mono", "IFN-gamma"), 0.05, thr)
        checks[f"DEGs with |log2FC| > {thr}"] = ((dd["direction"] != "ns").sum(), exp)
for name, (got, exp) in checks.items():
    ok = abs(int(got) - exp) <= 1
    print(f"  {'✓' if ok else '✗'} {name}: {int(got):,} (text: {exp:,})")
    if not ok:
        errors.append(("narrative numbers", name, f"{got} vs {exp}"))
P = L["PCA_ALL"]
for what, got, exp in (("PC1 %", P["var"][0], 66), ("PC2 %", P["var"][1], 17)):
    ok = abs(got - exp) < 1
    print(f"  {'✓' if ok else '✗'} global PCA {what}: {got:.1f} (text: {exp})")
    if not ok:
        errors.append(("narrative numbers", what, f"{got} vs {exp}"))

if "--fallback" in sys.argv:
    print("\nOn-demand DESeq2 (the cache is missing):")
    t = time.time()
    live = L["DEResults"](adata, None, verbose=False)
    r_live, r_cache = live.get("CD14_Mono", "IL-4"), L["DE"].get("CD14_Mono", "IL-4")
    both = r_live["tested"] & r_cache["tested"]
    r = np.corrcoef(r_live.loc[both, "log2fc"], r_cache.loc[both, "log2fc"])[0, 1]
    print(f"  {time.time() - t:.0f} s · {both.sum():,} genes · correlation of log2FC live vs cache: {r:.5f}")
    if r < 0.999:
        errors.append(("fallback", "DE", f"r = {r}"))

if "--execute" in sys.argv:
    import nbformat
    from nbclient import NotebookClient
    from build_interactive_notebooks import build
    print("\nExecuting the notebook from a fresh kernel …")
    nb = nbformat.read(build(LANG), 4)
    nb.cells.append(nbformat.v4.new_code_cell("assert not ERRORS, ERRORS\nprint('no error recorded while drawing the initial state of the exercises')"))
    t = time.time()
    NotebookClient(nb, timeout=1200, kernel_name="python3", resources={"metadata": {"path": str(ROOT)}}).execute()
    print(f"  all {sum(c.cell_type == 'code' for c in nb.cells)} code cells ran in {time.time() - t:.0f} s")

print(f"\n{n_calls} render calls, {state['figs']} figures drawn, figures in {FIGDIR}")
for name, vals, tb in errors[:10]:
    print(f"\n✗ {name}: {vals if not isinstance(vals, dict) else {k: v for k, v in vals.items() if k != 'samples'}}\n{tb}")
print(f"ERRORS: {len(errors)}")
sys.exit(1 if errors else 0)
