"""Generate every data figure and number used in slides/v3/main.tex.

The script *executes the code cells of CytokineDictionary_Practical_v3.ipynb* (so the slides always show exactly what the
practical computes: same objects, same numbers) and then draws each slide figure at its real printed size, with short titles,
legends placed in empty regions and fonts that stay readable on a 16:9 slide.

Usage (from anywhere):   python slides/v3/make_figures.py        (English: slides/v3)
                         python slides/v3/make_figures.py es     (Spanish: reads the _es notebook, writes slides/v3_es)
Needs:  pseudobulks_4ct_proteincoding.h5ad in the project root  (see the notebook for the download)
        pdftoppm (poppler) and inkscape, only to prepare the paper figure and the AnnData schema (skipped if already present)
Output: slides/v3/img/*.png  and  slides/v3/numbers.tex, bh_table.tex, libsize_table.tex
"""
import os, sys, json, subprocess, shutil, warnings
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LANG = sys.argv[1] if len(sys.argv) > 1 else "en"
def T(en, es):                                           # tiny translation helper for the text drawn inside the slide figures
    return es if LANG == "es" else en

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
HERE = HERE if LANG == "en" else HERE.parent / "v3_es"   # where the figures and tables of this language are written
OUT = HERE / "img"
OUT.mkdir(parents=True, exist_ok=True)
os.chdir(ROOT)                                           # the notebook uses a path relative to the project root
warnings.filterwarnings("ignore")

# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# 0. static assets: Figure 1 of the paper (panels a, b) and the AnnData schema
# ══════════════════════════════════════════════════════════════════════════════════════════════════════
def prepare_assets():
    paper = ROOT / "refs" / "Oesinghaus et al. - 2025 - A single-cell cytokine dictionary of human peripheral blood.pdf"
    target = OUT / "paper_fig1ab.png"
    if not target.exists() and paper.exists() and shutil.which("pdftoppm"):
        from PIL import Image
        tmp = OUT / "_paper_p5"
        subprocess.check_call(["pdftoppm", "-f", "5", "-l", "5", "-r", "300", "-png", str(paper), str(tmp)])
        page = Image.open(next(OUT.glob("_paper_p5*.png")))
        k = 300 / 110                                    # crop box measured on a 110-dpi render of the page
        crop = page.crop((int(105 * k), int(108 * k), int(522 * k), int(528 * k)))
        from PIL import ImageDraw
        dr = ImageDraw.Draw(crop)                          # hide the labels of panel c that peek in at the right edge
        W = crop.width
        for box in [(1100, 0, W, 70), (1052, 316, W, 358), (1052, 430, W, 470), (0, 1118, W, crop.height)]:
            dr.rectangle(box, fill="white")
        crop.save(target)
        for f in OUT.glob("_paper_p5*.png"): f.unlink()
    from PIL import Image
    def crop_page(pdf, page, dpi, box100, target_name):    # box measured on a 100-dpi render of the page
        target = OUT / target_name
        if target.exists() or not pdf.exists() or not shutil.which("pdftoppm"): return
        tmp = OUT / "_tmp"
        subprocess.check_call(["pdftoppm", "-f", str(page), "-l", str(page), "-r", str(dpi), "-png", str(pdf), str(tmp)])
        im = Image.open(next(OUT.glob("_tmp*.png"))); k = dpi / 100
        im.crop(tuple(int(v * k) for v in box100)).save(target)
        for f in OUT.glob("_tmp*.png"): f.unlink()
    refs = ROOT / "refs"
    tzec = next(refs.glob("Tzec*.pdf"), refs / "missing.pdf")
    crop_page(tzec, 3, 250, (185, 58, 650, 892), "tzec_fig1_bulk.png")                   # Fig. 1: bulk RNA-seq workflow
    crop_page(tzec, 4, 250, (183, 62, 645, 880), "tzec_fig2_sc.png")                     # Fig. 2: single-cell RNA-seq workflow
    crop_page(paper, 8, 250, (105, 98, 738, 368), "paper_fig2_donors.png")               # Fig. 2a-g of the cytokine paper: donor variability
    crop_page(paper, 5, 250, (95, 473, 492, 722), "paper_fig1de.png")                    # Fig. 1d-e: receptor expression and known cytokine targets
    _f = OUT / "paper_fig1de.png"
    if _f.exists():                                                                      # hide the labels of panel c that peek in at the right edge
        from PIL import ImageDraw
        _im = Image.open(_f).convert("RGB"); ImageDraw.Draw(_im).rectangle((int(0.925 * _im.width), 0, _im.width, int(0.60 * _im.height)), fill="white"); _im.save(_f)
    schema = ROOT / "refs" / "anndata_schema.svg"
    target = OUT / "anndata_schema.png"
    if not target.exists() and schema.exists() and shutil.which("inkscape"):
        subprocess.check_call(["inkscape", "--export-type=png", "--export-width=1600", f"--export-filename={target}", str(schema)],
                              stderr=subprocess.DEVNULL)
prepare_assets()
# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# 1. run the notebook code cells in one namespace (figures of the notebook are not needed here)
# ══════════════════════════════════════════════════════════════════════════════════════════════════════
nb = json.load(open(ROOT / ("CytokineDictionary_Practical_v3.ipynb" if LANG == "en" else "CytokineDictionary_Practical_v3_es.ipynb")))
ns = {"__name__": "__main__"}
plt.show = lambda *a, **k: plt.close("all")
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    src = "".join(c["source"])
    if "Install any missing dependencies" in src or not src.strip() or src.strip().startswith("# 🔧 your code here"):
        continue
    exec(compile(src, f"<cell {i}>", "exec"), ns)
    plt.close("all")
print("✅ notebook executed")
globals().update({k: v for k, v in ns.items() if not k.startswith("__")})

# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# 2. slide style: figures are drawn at ~8-10 inches wide and printed at ~14 cm (factor ~0.6)
# ══════════════════════════════════════════════════════════════════════════════════════════════════════
plt.rcParams.update({"font.size": 12, "axes.titlesize": 12.5, "axes.labelsize": 12, "xtick.labelsize": 11, "ytick.labelsize": 11,
                     "legend.fontsize": 10.5, "figure.dpi": 100, "savefig.dpi": 220, "axes.titleweight": "bold",
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
import matplotlib.ticker
rng = np.random.default_rng(7)

# every figure is drawn at FIG_SCALE × the size written in the code, so that fonts stay readable once the figure is
# shrunk to the width of a slide (the text sizes below are therefore "printed" sizes divided by ~0.6)
FIG_SCALE = 0.72
_figure_orig = plt.figure
def _scaled_figure(*a, **k):
    if k.get("figsize") is not None:
        w, h = k["figsize"]; k["figsize"] = (w * FIG_SCALE, h * FIG_SCALE)
    return _figure_orig(*a, **k)
plt.figure = _scaled_figure                                   # plt.subplots() calls plt.figure() internally

def save(fig, name):
    for ax_ in fig.axes:                                   # legends were excluded from tight_layout: include them in the saved area
        if ax_.get_legend() is not None: ax_.get_legend().set_in_layout(True)
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", facecolor="white", dpi=220)
    plt.close(fig); print("  saved", name)

obs = adata.obs
is_ctrl = (obs["treatment"] == "control").values
ct_arr = obs["cell_type"].astype(str).values

# ── Part 2 ──────────────────────────────────────────────────────────────────────────────────────────
# library size by cell type, and library size vs cells
FIG_SCALE = 0.66
fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6), gridspec_kw={"width_ratios": [1, 1.1]})
ax = axes[0]
for k, ct in enumerate(CELL_TYPES):
    m = ct_arr == ct
    v, vp = obs.loc[m & ~is_ctrl, "lib_size"].values, obs.loc[m & is_ctrl, "lib_size"].values
    ax.scatter(k - 0.12 + rng.uniform(-0.2, 0.2, len(v)), v, s=6, alpha=0.4, linewidths=0, color=CT_COLORS[ct], marker=CT_MARKERS[ct])
    ax.scatter(k + 0.34 + rng.uniform(-0.06, 0.06, len(vp)), vp, s=26, facecolors="none", edgecolors="k", linewidths=0.9)
    ax.hlines(np.median(v), k - 0.42, k + 0.18, color="k", lw=2)
ax.set_yscale("log"); ax.set_ylim(5e4, 5e8)
ax.set_xticks(range(4)); ax.set_xticklabels(["CD4", "CD8", "NK", "Mono"]); ax.set_ylabel(T("library size", "tamaño de librería"))
ax.set_title(T("by cell type (○ = PBS)", "por tipo celular (○ = PBS)"))
ax = axes[1]
for ct in CELL_TYPES:
    m = (ct_arr == ct) & ~is_ctrl
    ax.scatter(obs.loc[m, "psbulk_n_cells"], obs.loc[m, "lib_size"], s=6, alpha=0.4, linewidths=0, color=CT_COLORS[ct], marker=CT_MARKERS[ct], label=CT_SHORT[ct])
ax.scatter(obs.loc[is_ctrl, "psbulk_n_cells"], obs.loc[is_ctrl, "lib_size"], s=26, facecolors="none", edgecolors="k", linewidths=0.9, label="PBS")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(T("cells in the pseudobulk", "células en el pseudobulk")); ax.set_ylabel(T("library size", "tamaño de librería"))
ax.set_title(T("follows the number of cells", "sigue al nº de células"))
plt.tight_layout(); save(fig, "p2_libsize")

FIG_SCALE = 0.64
# ACTB before/after CPM
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.3), sharex=True)
m_mono = (obs["cell_type"] == "CD14_Mono").values
j = adata.var_names.get_loc("ACTB")                                           # position in the FILTERED object: use its own layers
libsz = adata.obs["lib_size"].values
for ax, vals, label in zip(axes, [adata.layers["counts"][m_mono, j], adata.layers["cpm"][m_mono, j]], [T("raw count", "conteo crudo"), "CPM"]):
    ax.scatter(libsz[m_mono], vals, s=7, alpha=0.5, color=CT_COLORS["CD14_Mono"], linewidths=0)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    r = np.corrcoef(np.log10(libsz[m_mono]), np.log10(vals + 1))[0, 1]
    ax.set_title(f"ACTB: {label}\nr = {r:+.2f}"); ax.set_xlabel(T("library size of the sample", "tamaño de librería")); ax.set_ylabel(label)
plt.tight_layout(); save(fig, "p2_actb")

FIG_SCALE = 0.58
# CPM vs log: variance against mean
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.4))
ok = (mean_cpm > 0.01) & (var_cpm > 0)
axes[0].scatter(mean_cpm[ok], var_cpm[ok], s=3, alpha=0.25, color="#2a78d6", linewidths=0, rasterized=True)
xx = np.array([0.01, 1e5]); axes[0].plot(xx, xx ** 2 * 0.5, color="k", ls="--", lw=1.3, label=T("slope 2", "pendiente 2"))
axes[0].set_xscale("log"); axes[0].set_yscale("log"); axes[0].legend(loc="lower right")
axes[0].set_xlabel(T("mean CPM", "CPM medio")); axes[0].set_ylabel(T("variance", "varianza")); axes[0].set_title("CPM")
ok = var_log > 0
axes[1].scatter(mean_log[ok], var_log[ok], s=3, alpha=0.25, color="#2a78d6", linewidths=0, rasterized=True)
axes[1].set_xlabel(T("mean log₂(CPM + 1)", "log₂(CPM + 1) medio")); axes[1].set_ylabel(T("variance", "varianza")); axes[1].set_title("log₂(CPM + 1)")
plt.tight_layout(); save(fig, "p2_logvar")

FIG_SCALE = 0.66
# the toy example
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.5))
genes_t = ["A", "B", "C", "D"]
raw_t = np.array([[1000, 200], [50, 40], [2500, 500], [1450, 260]])
cpm_t = raw_t / raw_t.sum(axis=0) * 1e6 / 1e3                                  # in thousands of CPM
for ax, M, title, ylab in [(axes[0], raw_t, T("raw counts", "conteos crudos"), T("counts", "conteos")), (axes[1], cpm_t, "CPM", T("CPM (thousands)", "CPM (miles)"))]:
    x = np.arange(4); w = 0.38
    ax.bar(x - w / 2, M[:, 0], w, color=PBS_COLOR, label=T("PBS (5,000 counts in total)", "PBS (5.000 conteos en total)"))
    ax.bar(x + w / 2, M[:, 1], w, color="#eb6834", label=T("IFN-γ (1,000 counts in total)", "IFN-γ (1.000 conteos en total)"))
    for xi in range(4):
        ratio = M[xi, 1] / M[xi, 0]
        ax.text(xi, max(M[xi]) + 0.02 * M.max(), f"×{ratio:.1f}", ha="center", va="bottom", fontsize=12, fontweight="bold", color=UP_COLOR if ratio > 1.05 else "#52514e")
    ax.set_xticks(x); ax.set_xticklabels(genes_t); ax.set_ylabel(ylab); ax.set_title(title, pad=8); ax.set_ylim(0, M.max() * 1.22)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=[0, 0.09, 1, 1]); save(fig, "p2_toy")

FIG_SCALE = 0.68
# Poisson noise: the same gene measured twice
fig, axes = plt.subplots(1, 4, figsize=(9.6, 2.7), sharex=True, sharey=True)
for ax, lam in zip(axes, [2, 20, 200, 2000]):
    a, b = rng.poisson(lam, 5000), rng.poisson(lam, 5000)
    lfc = np.log2((a + 1) / (b + 1))
    ax.hist(lfc, bins=np.linspace(-4, 4, 41), color="#2a78d6" if lam >= 200 else PBS_COLOR)
    ax.set_title(T(f"mean {lam}\n{np.mean(np.abs(lfc) > 1):.0%} off by >2×", f"media {lam}\n>2×: {np.mean(np.abs(lfc) > 1):.0%}"), fontsize=10)
    ax.set_xlabel("log₂FC"); ax.set_xticks([-3, 0, 3])
axes[0].set_ylabel("genes")
plt.tight_layout(); save(fig, "p2_poisson")

FIG_SCALE = 0.62
# the filter
n_ok = (counts >= MIN_COUNT).sum(axis=0)                                      # `counts`: raw counts of the 18,120 genes (before filtering)
fig, ax = plt.subplots(figsize=(6.4, 2.9))
ax.hist(np.log10(n_ok + 1), bins=60, color="#2a78d6")
ax.axvline(np.log10(MIN_SAMPLES + 1), color=UP_COLOR, ls="--", lw=1.8)
ax.annotate(T("12 samples", "12 muestras"), xy=(np.log10(MIN_SAMPLES + 1), ax.get_ylim()[1] * 0.85), xytext=(7, 0), textcoords="offset points", ha="left", color=UP_COLOR, fontsize=11.5, fontweight="bold")
ax.set_xlabel(T("log₁₀(1 + number of samples with ≥ 10 counts)", "log₁₀(1 + nº de muestras con ≥ 10 conteos)")); ax.set_ylabel("genes")
plt.tight_layout(); save(fig, "p2_filter")

# block-mean heatmaps (cell type × condition, mean over the 12 donors)
FIG_SCALE = 0.68
def slide_heatmap(gene_sets, name, conditions=CONDITIONS):
    genes = [g for gs in gene_sets.values() for g in gs]
    sizes = [len(gs) for gs in gene_sets.values()]
    cols  = [(ct, cy) for ct in CELL_TYPES for cy in conditions]
    def blockmean(layer):
        out = np.zeros((len(genes), len(cols)))
        for jj, (ct, cy) in enumerate(cols):
            sub = adata[[f"{d}|{ct}|{cy}" for d in DONORS], genes]
            out[:, jj] = np.asarray(sub.layers["counts"] if layer == "counts" else sub.X).mean(axis=0)
        return out
    libb = np.array([np.log10(adata.obs.loc[[f"{d}|{ct}|{cy}" for d in DONORS], "lib_size"].values).mean() for ct, cy in cols])[None, :]
    fig = plt.figure(figsize=(9.8, 0.3 * len(genes) + 2.5))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.8, 0.8, len(genes)], width_ratios=[1, 1, 0.03], hspace=0.06, wspace=0.12, top=0.88, bottom=0.2, left=0.1, right=0.93)
    for c, (title, layer) in enumerate([(T("BEFORE: raw counts", "ANTES: crudos"), "counts"), (T("AFTER: log₂(CPM + 1)", "DESPUÉS: log₂(CPM+1)"), "x")]):
        M = blockmean(layer); Z = (M - M.mean(1, keepdims=True)) / M.std(1, ddof=1, keepdims=True)
        ax_ct, ax_lib, ax = (fig.add_subplot(gs[r, c]) for r in range(3))
        for k, ct in enumerate(CELL_TYPES):
            ax_ct.add_patch(plt.Rectangle((k * len(conditions) - 0.5, 0), len(conditions), 1, color=CT_COLORS[ct], lw=0))
            ax_ct.text(k * len(conditions) + (len(conditions) - 1) / 2, 0.5, CT_SHORT[ct], ha="center", va="center", fontsize=10.5, fontweight="bold",
                       color="black" if ct == "NK_CD56low" else "white")
        ax_ct.set_xlim(-0.5, len(cols) - 0.5); ax_ct.set_ylim(0, 1); ax_ct.axis("off"); ax_ct.set_title(title, fontsize=14, pad=8)
        ax_lib.imshow(libb, aspect="auto", cmap="Greys", interpolation="nearest")
        ax_lib.set_xticks([]); ax_lib.set_yticks([0] if c == 0 else [])
        if c == 0: ax_lib.set_yticklabels([T("library\nsize", "tamaño\nlibrería")], fontsize=10.5)
        for sp in ax_lib.spines.values(): sp.set_visible(False)
        im = ax.imshow(Z, aspect="auto", cmap=DIVERGING, vmin=-2, vmax=2, interpolation="nearest")
        ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes if c == 0 else [], fontsize=12)
        ax.set_xticks(range(len(cols))); ax.set_xticklabels([cy.replace("-gamma", "-γ").replace("-beta", "-β") for _, cy in cols], rotation=90, fontsize=10)
        for b in np.cumsum(sizes)[:-1] - 0.5: ax.axhline(b, color="white", lw=2)
        for k in range(1, len(CELL_TYPES)): ax.axvline(k * len(conditions) - 0.5, color="white", lw=2)
        for sp in ax.spines.values(): sp.set_visible(False)
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[2, 2])); cb.set_label("z-score", fontsize=12); cb.ax.tick_params(labelsize=10.5)
    save(fig, name)

slide_heatmap(CELL_TYPE_MARKERS, "p2_heat_markers")
slide_heatmap(IMMUNE_GENES, "p2_heat_immune")

# ── Part 3 ──────────────────────────────────────────────────────────────────────────────────────────
# two genes: rotation and projection  (P, Pc, V, S, mu, is_ifn come from the notebook's section 3.1)
FIG_SCALE = 0.72
scores2g, var2g = Pc @ V, S ** 2 / np.sum(S ** 2) * 100
fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.9), gridspec_kw={"width_ratios": [1, 1.25]})
ax = axes[0]
for flag, col, lab in [(False, NS_COLOR, T("other", "otras")), (True, UP_COLOR, T("interferons", "interferones"))]:
    ax.scatter(P[is_ifn == flag, 0], P[is_ifn == flag, 1], s=9, alpha=0.65, color=col, linewidths=0, label=lab)
for i in rng.choice(len(P), 12, replace=False):
    proj = mu + (Pc[i] @ V[:, 0]) * V[:, 0]
    ax.plot([P[i, 0], proj[0]], [P[i, 1], proj[1]], color="k", lw=0.9, ls=":"); ax.scatter(*proj, s=14, color="k", zorder=5)
for k, nm in enumerate(["PC1", "PC2"]):
    d_ = V[:, k] * 1.5 * S[k] / np.sqrt(len(P) - 1)
    ax.annotate("", xy=mu + d_, xytext=mu - d_, arrowprops=dict(arrowstyle="<|-|>", color="k", lw=1.8), annotation_clip=False)
    ax.annotate(nm, xy=mu + d_, xytext=(6, 4), textcoords="offset points", fontsize=12, fontweight="bold", annotation_clip=False)
ax.set_aspect("equal"); ax.set_xlabel("ISG15  (log₂ CPM + 1)"); ax.set_ylabel("MX1  (log₂ CPM + 1)"); ax.legend(loc="lower right", markerscale=2)
ax.set_title(T("original axes", "ejes originales"))
ax = axes[1]
for flag, col in [(False, NS_COLOR), (True, UP_COLOR)]:
    ax.scatter(scores2g[is_ifn == flag, 0], scores2g[is_ifn == flag, 1], s=9, alpha=0.65, color=col, linewidths=0)
ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6); ax.set_aspect("equal")
ax.set_xlabel(T(f"PC1 ({var2g[0]:.0f} % of the variance)", f"PC1 ({var2g[0]:.0f} % de la varianza)")); ax.set_ylabel(f"PC2 ({var2g[1]:.0f} %)"); ax.set_title(T("rotated axes", "ejes rotados"))
plt.tight_layout(); save(fig, "p3_twogenes")

# the PCA of all samples, shown in two figures of two panels each
def pca_panel(ax, data, color_by, pcs=(0, 1), s=10, legend_kw=None, title=None):
    sc, (a, b), var, ob = data.obsm["X_pca"], pcs, data.uns["pca_var"], data.obs
    kw = dict(linewidths=0.25, edgecolors="white")
    if color_by == "cell_type":
        for ct in CELL_TYPES:
            m = (ob["cell_type"] == ct).values
            ax.scatter(sc[m, a], sc[m, b], s=s, alpha=0.8, color=CT_COLORS[ct], marker=CT_MARKERS[ct], label=CT_LABEL[ct], **kw)
    elif color_by == "donor":
        for i, d_ in enumerate(DONORS):
            m = (ob["donor"] == d_).values
            ax.scatter(sc[m, a], sc[m, b], s=s, alpha=0.85, color=DONOR_HUES[i % 6], marker="o" if i < 6 else "^", label=d_.replace("Donor", "D"), **kw)
    elif color_by == "cytokine":
        m = ~ob["cytokine"].isin(["PBS"] + CYTOKINES_PANEL).values
        ax.scatter(sc[m, a], sc[m, b], s=s * 0.6, color=NS_COLOR, alpha=0.5, linewidths=0, label=T("others", "otras"), zorder=1)
        for i, cy in enumerate(CYTOKINES_PANEL):
            m = (ob["cytokine"] == cy).values
            ax.scatter(sc[m, a], sc[m, b], s=s * 1.5, color=CYT_HUES[i], marker=CYT_MARKERS[i], label=cy, zorder=2, **kw)
        m = (ob["cytokine"] == "PBS").values
        ax.scatter(sc[m, a], sc[m, b], s=s * 2.6, color="black", label="PBS", zorder=3, **kw)
    elif color_by == "treatment":
        m = (ob["treatment"] == "treated").values
        ax.scatter(sc[m, a], sc[m, b], s=s * 0.8, color=NS_COLOR, alpha=0.6, linewidths=0, label=T(f"treated ({m.sum():,})", f"tratadas ({m.sum():,})"), zorder=1)
        ax.scatter(sc[~m, a], sc[~m, b], s=s * 2.6, color="black", label=T(f"control PBS ({(~m).sum()})", f"control PBS ({(~m).sum()})"), zorder=2, **kw)
    if legend_kw is not None:                              # legends go BELOW the plot: they can never hide the data
        leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), borderaxespad=4.4, markerscale=1.4, handletextpad=0.2, labelspacing=0.25,
                        columnspacing=0.9, fontsize=9.5, **legend_kw)
        leg.set_in_layout(False)                           # tight_layout ignores it (bbox_inches="tight" still saves it)
    ax.set_xlabel(f"PC{a + 1} ({var[a]:.1f} %)"); ax.set_ylabel(f"PC{b + 1} ({var[b]:.1f} %)")
    ax.set_title(title or color_by)

FIG_SCALE = 0.68
fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.5))
pca_panel(axes[0], adata, "cell_type", title=T("by cell type", "por tipo celular"), legend_kw=dict(ncol=2))
pca_panel(axes[1], adata, "donor", title=T("by donor", "por donante"), legend_kw=dict(ncol=4))
plt.tight_layout(); save(fig, "p3_all_ct_donor")

FIG_SCALE = 0.68
fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.5))
pca_panel(axes[0], adata, "cytokine", title=T("by cytokine", "por citoquina"), legend_kw=dict(ncol=4))
pca_panel(axes[1], adata, "treatment", title=T("by control vs treated", "control vs tratadas"), legend_kw=dict(ncol=1))
plt.tight_layout(); save(fig, "p3_all_cyt_trt")

FIG_SCALE = 0.66
# scree plot + R² table (side by side)
R2 = pc_r2_table(adata, 8); var8 = adata.uns["pca_var"][:8]
fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.6), gridspec_kw={"width_ratios": [0.9, 1.5]})
axes[0].barh(range(8), var8, color="#2a78d6"); axes[0].invert_yaxis(); axes[0].set_yticks(range(8)); axes[0].set_yticklabels(R2.index)
for i, v in enumerate(var8): axes[0].text(v + 1, i, f"{v:.1f}", va="center", fontsize=10)
axes[0].set_xlabel(T("variance explained (%)", "varianza explicada (%)")); axes[0].set_title(T("how much?", "¿cuánta?"))
im = axes[1].imshow(R2.values, cmap=SEQUENTIAL, vmin=0, vmax=1, aspect="auto")
for i in range(R2.shape[0]):
    for jj in range(R2.shape[1]):
        v = R2.values[i, jj]; axes[1].text(jj, i, f"{v:.2f}", ha="center", va="center", fontsize=10, color="w" if v > 0.55 else "k")
axes[1].set_xticks(range(R2.shape[1])); axes[1].set_xticklabels([T("cell\ntype", "tipo\ncelular"), T("donor", "donante"), T("cyto-\nkine", "cito-\nquina"), T("ctrl vs\ntreated", "ctrl vs\ntratado"), T("library\nsize", "tamaño\nlibrería")], fontsize=9)
axes[1].set_yticks(range(8)); axes[1].set_yticklabels(R2.index); axes[1].set_title(T("what? (R² with each variable)", "¿qué? (R² con cada variable)"))
plt.colorbar(im, ax=axes[1], shrink=0.85, pad=0.02)
plt.tight_layout(); save(fig, "p3_r2")

# loadings of PC1-PC3
def loading_axes(axs, load, genes, pcs, var, n=6):
    for ax, k in zip(np.atleast_1d(axs), pcs):
        s_ = pd.Series(load[:, k], index=genes).sort_values(); top = pd.concat([s_.head(n), s_.tail(n)])
        ax.barh(top.index, top.values, color=[DOWN_COLOR if v < 0 else UP_COLOR for v in top.values], height=0.8)
        ax.axvline(0, color="k", lw=0.6); ax.tick_params(axis="y", labelsize=9.5); ax.xaxis.set_major_locator(plt.MaxNLocator(2))
        ax.set_xlabel(T("loading", "carga")); ax.set_title(f"PC{k + 1} ({var[k]:.1f} %)")
FIG_SCALE = 0.68
fig, axes = plt.subplots(1, 3, figsize=(9.8, 3.7))
loading_axes(axes, loadings, adata.var_names[hvg], (0, 1, 2), adata.uns["pca_var"], n=5)
plt.tight_layout(); save(fig, "p3_loadings")

FIG_SCALE = 0.68
# monocytes only
fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.5))
pca_panel(axes[0], mono, "donor", title=T("monocytes: by donor", "monocitos: por donante"), legend_kw=dict(ncol=3))
pca_panel(axes[1], mono, "cytokine", title=T("monocytes: by cytokine", "monocitos: por citoquina"), legend_kw=dict(ncol=2))
plt.tight_layout(); save(fig, "p3_mono_donor_cyt")

FIG_SCALE = 0.68
fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.7), gridspec_kw={"width_ratios": [1.15, 1]})
pca_panel(axes[0], mono, "cytokine", pcs=(2, 3), title=T("monocytes: PC3 vs PC4", "monocitos: PC3 vs PC4"), legend_kw=dict(ncol=3))
loading_axes(axes[1], load_m, mono.var_names[hvg_m], (3,), mono.uns["pca_var"], n=4)
plt.tight_layout(); save(fig, "p3_mono_pc34")

# ── Part 4 ──────────────────────────────────────────────────────────────────────────────────────────
FIG_SCALE = 0.62
fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5))
paired_plot(axes[0], x_pbs, x_cyt, EXAMPLE_GENE, CYTOKINE)
paired_plot(axes[1], ctrl_pbs, ctrl_cyt, CONTROL_GENE, CYTOKINE)
for ax in axes: ax.title.set_fontsize(10.5)
plt.tight_layout(); save(fig, "p4_paired")

FIG_SCALE = 0.6
xx = np.linspace(-8, 21, 600)
fig, ax = plt.subplots(figsize=(5.8, 3.1))
ax.plot(xx, stats.t.pdf(xx, df=11), color="k", lw=1.6)
for sgn in (1, -1):
    tail = xx[sgn * xx >= abs(t_obs)]; ax.fill_between(tail, stats.t.pdf(tail, df=11), color=UP_COLOR, alpha=0.7)
ax.axvline(t_obs, color=UP_COLOR, ls="--")
ax.annotate(f"t = {t_obs:.1f}\np = {p_val:.0e}", xy=(t_obs, 0.16), xytext=(8, 0), textcoords="offset points", ha="left", color=UP_COLOR, fontsize=10, fontweight="bold")
ax.set_xlim(-8, 21); ax.set_xlabel(T("t statistic", "estadístico t")); ax.set_ylabel(T("density", "densidad")); ax.set_title(T("t distribution, 11 degrees of freedom", "distribución t, 11 grados de libertad"))
plt.tight_layout(); save(fig, "p4_tdist")

# DESeq2 under the hood: size factors and dispersions
FIG_SCALE = 0.66
sf_ = pd.Series(de.attrs["size_factors"]); lib24 = adata.obs.loc[sf_.index, "lib_size"].values
is_ref_ = np.array([n.endswith("|PBS") for n in sf_.index])
fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6))
ax = axes[0]
ax.scatter(lib24[is_ref_], sf_.values[is_ref_], s=40, color=PBS_COLOR, label="PBS"); ax.scatter(lib24[~is_ref_], sf_.values[~is_ref_], s=40, color=CYTOKINE_COLOR, label=CYTOKINE)
ax.set_xscale("log"); ax.set_yscale("log"); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_yticks([0.2, 0.5, 1, 2, 5]); ax.set_yticklabels(["0.2", "0.5", "1", "2", "5"])
ax.set_xlabel(T("library size", "tamaño de librería")); ax.set_ylabel(T("size factor", "factor de tamaño")); ax.set_title(T("size factors", "factores de tamaño")); ax.legend(loc="upper left")
ax = axes[1]
dd = de[de["tested"] & (de["disp_genewise"] > 1e-6)]
ax.scatter(dd["base_mean"], dd["disp_genewise"], s=3, alpha=0.4, color=PBS_COLOR, linewidths=0, rasterized=True)
ax.scatter(dd["base_mean"], dd["disp_final"], s=3, alpha=0.25, color="#2a78d6", linewidths=0, rasterized=True)
o_ = np.argsort(dd["base_mean"].values); ax.plot(dd["base_mean"].values[o_], dd["disp_trend"].values[o_], color=UP_COLOR, lw=2)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(T("mean of normalized counts", "media de conteos normalizados")); ax.set_ylabel(T("dispersion (noise)", "dispersión (ruido)")); ax.set_title(T("dispersion", "dispersión"))
leg = ax.legend(handles=[Line2D([], [], marker="o", color=PBS_COLOR, ls="", label=T("gene-wise", "por gen")), Line2D([], [], marker="o", color="#2a78d6", ls="", label=T("shrunk", "contraída")),
                   Line2D([], [], color=UP_COLOR, lw=2, label=T("trend", "tendencia"))], loc="upper center", bbox_to_anchor=(0.5, 0.0), borderaxespad=4.4, ncol=3, fontsize=9)
leg.set_in_layout(False)
plt.tight_layout(); save(fig, "p4_deseq_diag")

# DESeq2 against the log-CPM paired t-test
FIG_SCALE = 0.66
fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6))
ax = axes[0]
ax.scatter(de_t.loc[both, "log2fc"], de.loc[both, "log2fc"], s=3, alpha=0.3, color="#2a78d6", linewidths=0, rasterized=True)
lim = 4.5; ax.plot([-lim, lim], [-lim, lim], color="k", lw=0.8, ls="--"); ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
r_eff = np.corrcoef(de_t.loc[both, "log2fc"], de.loc[both, "log2fc"])[0, 1]
ax.set_xlabel(T("log2FC, log-CPM t-test", "log2FC, prueba t log-CPM")); ax.set_ylabel("log2FC, DESeq2"); ax.set_title(T(f"effect sizes (r = {r_eff:.2f})", f"tamaños de efecto (r = {r_eff:.2f})"))
ax = axes[1]
bins_ = np.linspace(-3, 3, 90)
for res_, col, lab in [(de_t, PBS_COLOR, T("log-CPM t-test", "prueba t log-CPM")), (de, "#2a78d6", "DESeq2")]:
    v_ = res_.loc[res_["tested"], "log2fc"]; ax.hist(v_, bins=bins_, histtype="stepfilled", alpha=0.55, color=col, label=f"{lab} ({T('median', 'mediana')} {v_.median():+.2f})")
ax.axvline(0, color="k", lw=0.8); ax.set_xlabel("log₂ fold change"); ax.set_ylabel("genes"); ax.set_title(T("distribution of the log₂FC", "distribución del log₂FC"))
leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), borderaxespad=4.4, fontsize=9); leg.set_in_layout(False)
plt.tight_layout(); save(fig, "p4_vs_ttest")

FIG_SCALE = 0.6
fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.3), sharey=True)
p_real_ = de.loc[de["tested"], "pval"].values
for ax, p_, title, col in [(axes[0], p_real_, T("real data", "datos reales"), "#2a78d6"), (axes[1], p_null, T("null (no effect)", "nula (sin efecto)"), PBS_COLOR)]:
    ax.hist(p_, bins=40, range=(0, 1), color=col); ax.axhline(len(p_) / 40, color="k", ls="--", lw=0.9)
    ax.set_xlabel(T("p-value", "valor p")); ax.set_xticks([0, 0.5, 1]); ax.set_title(T(f"{title}\n{(p_ < 0.05).sum():,} with p < 0.05", f"{title}\n{(p_ < 0.05).sum():,} con p < 0,05".replace(",05", ".05")), fontsize=10.5)
axes[0].set_ylabel("genes")
plt.tight_layout(); save(fig, "p4_pvalues")

# volcano plots (labels at a slide-friendly size)
def slide_volcano(ax, dfr, title, n_label=8, fs=9.5, xlabel=None):
    d_ = dfr[dfr["tested"]].copy(); d_["nlp"] = -np.log10(d_["padj"].clip(lower=1e-300))
    for k, col, size, alpha in [("ns", NS_COLOR, 6, 0.5), ("down", DOWN_COLOR, 12, 0.8), ("up", UP_COLOR, 12, 0.8)]:
        sub = d_[d_["direction"] == k]; ax.scatter(sub["log2fc"], sub["nlp"], s=size, color=col, alpha=alpha, linewidths=0, rasterized=True)
    for x in (-LFC_THRESH, LFC_THRESH): ax.axvline(x, color="#999", lw=0.8, ls="--")
    ax.axhline(-np.log10(PADJ_THRESH), color="#999", lw=0.8, ls="--")
    ax.margins(y=0.12); x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    sig = d_[d_["direction"] != "ns"].copy(); sig["score"] = sig["nlp"] * sig["log2fc"].abs()
    sig = sig[~((sig["nlp"] > y0 + 0.85 * (y1 - y0)) & ((sig["log2fc"] < x0 + 0.3 * (x1 - x0)) | (sig["log2fc"] > x1 - 0.3 * (x1 - x0))))]
    top = pd.concat([sig[sig["direction"] == "up"].nlargest(n_label, "score"), sig[sig["direction"] == "down"].nlargest(n_label, "score")]) if n_label else sig.iloc[:0]
    texts = [ax.text(r.log2fc, r.nlp, g, fontsize=fs) for g, r in top.iterrows()]
    if HAS_ADJUSTTEXT and texts: adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#999", lw=0.5))
    ax.annotate(f"↑ {(d_['direction'] == 'up').sum():,}", xy=(0.98, 0.98), xycoords="axes fraction", ha="right", va="top", color=UP_COLOR, fontsize=13, fontweight="bold")
    ax.annotate(f"↓ {(d_['direction'] == 'down').sum():,}", xy=(0.02, 0.98), xycoords="axes fraction", ha="left", va="top", color=DOWN_COLOR, fontsize=13, fontweight="bold")
    ax.set_xlabel(xlabel or T("log₂ fold change", "log₂ fold change")); ax.set_ylabel(T("−log₁₀ adj. p-value", "−log₁₀ p ajustado")); ax.set_title(title)

FIG_SCALE = 0.7
fig, ax = plt.subplots(figsize=(6.0, 4.5))
slide_volcano(ax, de, f"{CYTOKINE} vs PBS · {T('monocytes', 'monocitos')}", n_label=4, fs=9)
plt.tight_layout(); save(fig, "p4_volcano")

# two cytokines in the first cell type of the comparison: shared / specific responses and the direct comparison
FIG_SCALE = 0.58
ct0 = CELL_TYPES_COMPARE[0]; A0, B0, AB0 = cmp[ct0]["A"], cmp[ct0]["B"], cmp[ct0]["AB"]
ok0 = A0["tested"] & B0["tested"]; sa, sb = (A0["direction"] != "ns") & ok0, (B0["direction"] != "ns") & ok0
fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))
ax = axes[0]
for mk, col, lab in [(ok0 & ~sa & ~sb, NS_COLOR, T("neither", "ninguna")), (ok0 & sa & ~sb, "#2a78d6", T(f"{CYTOKINE_A} only", f"solo {CYTOKINE_A}")),
                     (ok0 & ~sa & sb, "#eb6834", T(f"{CYTOKINE_B} only", f"solo {CYTOKINE_B}")), (ok0 & sa & sb, "#6250d6", T("both", "ambas"))]:
    ax.scatter(A0.loc[mk, "log2fc"], B0.loc[mk, "log2fc"], s=4, alpha=0.55, color=col, linewidths=0, label=f"{lab} ({mk.sum():,})", rasterized=True)
ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel(f"log2FC {CYTOKINE_A} vs PBS"); ax.set_ylabel(f"log2FC {CYTOKINE_B} vs PBS"); ax.set_title(T("shared and specific responses", "compartidas y específicas"))
leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), borderaxespad=4.4, markerscale=3, fontsize=8.5, ncol=2); leg.set_in_layout(False)
slide_volcano(axes[1], AB0, f"{CYTOKINE_A} vs {CYTOKINE_B}", n_label=3, fs=9, xlabel=f"log₂FC ({CYTOKINE_A} / {CYTOKINE_B})")
plt.tight_layout(); save(fig, "p4_compare")

# Benjamini–Hochberg: sorted p-values against the line alpha·k/m
ps = np.sort(de.loc[de["tested"], "pval"].values); mm = len(ps); kk = np.arange(1, mm + 1)
kstar = int(np.where(ps <= 0.05 * kk / mm)[0].max() + 1)
FIG_SCALE = 0.62
fig, ax = plt.subplots(figsize=(5.6, 3.5))
ax.plot(kk, ps, color="#2a78d6", lw=2.5, label=T("p-values", "valores p"))
ax.plot(kk, 0.05 * kk / mm, color=UP_COLOR, lw=2, ls="--", label=r"line $\alpha\,k/m$".replace("line", T("line", "recta")))
ax.axvline(kstar, color="k", lw=0.8)
ax.annotate(f"k* = {kstar:,}", xy=(kstar, 0.006), xytext=(8, 0), textcoords="offset points", ha="left", va="bottom", fontsize=10.5, fontweight="bold")
ax.set_xlim(0, mm); ax.set_ylim(0, 0.06); ax.set_xlabel(T("rank k", "rango k")); ax.set_ylabel(T("p-value", "valor p")); ax.legend(loc="upper left", fontsize=9, handlelength=1.3)
plt.tight_layout(); save(fig, "p4_bh")

# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# 3. numbers and tables used verbatim in the slides
# ══════════════════════════════════════════════════════════════════════════════════════════════════════
def tex_int(v): return f"{int(round(v)):,}".replace(",", "\\,")
treated = obs["treatment"] == "treated"
med = obs[treated].groupby("cell_type", observed=True)[["psbulk_n_cells", "umis_per_cell", "lib_size"]].median()
sig_d_ = (de["padj"] < PADJ_THRESH) & (de["log2fc"].abs() > LFC_THRESH) & de["tested"]
sig_t_ = (de_t["padj"] < PADJ_THRESH) & (de_t["log2fc"].abs() > LFC_THRESH) & de_t["tested"]
nums = {
    "nsamples": tex_int(4224), "ngenesraw": tex_int(18120), "ngenes": tex_int(adata.n_vars), "ngenesremoved": tex_int(18120 - adata.n_vars),
    "umiMono": tex_int(med.loc["CD14_Mono", "umis_per_cell"]), "umiNK": tex_int(med.loc["NK_CD56low", "umis_per_cell"]),
    "pbsratio": f"{obs.loc[~treated, 'lib_size'].median() / obs.loc[treated, 'lib_size'].median():.1f}",
    "toptwentycpm": f"{top20_cpm:.0%}".replace("%", "\\%"), "toptwentylog": f"{top20_log:.1%}".replace("%", "\\%"),
    "pcaone": f"{adata.uns['pca_var'][0]:.0f}", "pcatwo": f"{adata.uns['pca_var'][1]:.0f}", "pcathree": f"{adata.uns['pca_var'][2]:.1f}",
    "pcafirstfour": f"{adata.uns['pca_var'][:4].sum():.0f}",
    "ntested": tex_int(de["tested"].sum()), "nup": tex_int((de["direction"] == "up").sum()), "ndown": tex_int((de["direction"] == "down").sum()),
    "nrawsig": tex_int((de.loc[de["tested"], "pval"] < 0.05).sum()), "nadjsig": tex_int((de.loc[de["tested"], "padj"] < 0.05).sum()),
    "nexpectedfp": tex_int(0.05 * de["tested"].sum()), "nnullraw": tex_int((p_null < 0.05).sum()), "nnulladj": tex_int((padj_null < 0.05).sum()),
    "tobs": f"{t_obs:.1f}", "pobs": f"{p_val:.1e}".replace("e-0", "\\times10^{-").replace("e-", "\\times10^{-") + "}",
    "lfcobs": f"{mean_d:.2f}", "foldobs": f"{2 ** mean_d:.1f}", "kstar": tex_int(kstar),
    "mediancpm": "$" + f"{de_t.loc[de_t['tested'], 'log2fc'].median():+.2f}" + "$", "mediandeseq": "$" + f"{de.loc[de['tested'], 'log2fc'].median():+.2f}" + "$",
    "ndegdeseq": tex_int(sig_d_.sum()), "ndegt": tex_int(sig_t_.sum()), "ndegboth": tex_int((sig_d_ & sig_t_).sum()), "reffect": f"{r_eff:.2f}",
}
with open(HERE / "numbers.tex", "w") as f:
    for k, v in nums.items(): f.write(f"\\newcommand{{\\num{k}}}{{{v}}}\n")

with open(HERE / "bh_table.tex", "w") as f:
    f.write("\\begin{tabular}{lrrrr}\\toprule\n$i$ & $p_{(i)}$ & $p_{(i)}\\,m/i$ & $p^{\\mathrm{adj}}_{(i)}$ & " + T("sig.", "sig.") + "\\\\\\midrule\n")
    for i in range(10):
        f.write(f"{i + 1} & {p_toy[i]:.3f} & {step1[i]:.3f} & {bh[i]:.3f} & {'\\checkmark' if bh[i] < 0.05 else ''}\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n")

with open(HERE / "libsize_table.tex", "w") as f:
    f.write("\\begin{tabular}{lrrr}\\toprule\n & " + T("cells", "células") + " & " + T("UMIs / cell", "UMI / célula") + " & " + T("library size", "librería") + "\\\\\\midrule\n")
    for ct, lab in [("CD4_T_cell", "CD4 T"), ("CD8_T_cell", "CD8 T"), ("NK_CD56low", "NK"), ("CD14_Mono", T("Monocytes", "Monocitos"))]:
        f.write(f"{lab} & {tex_int(med.loc[ct, 'psbulk_n_cells'])} & {tex_int(med.loc[ct, 'umis_per_cell'])} & {med.loc[ct, 'lib_size'] / 1e6:.1f}\\,M\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n")

# the two-cytokine comparison, for every cell type of CELL_TYPES_COMPARE
with open(HERE / "degcount_table.tex", "w") as f:
    f.write("\\begin{tabular}{@{}lrrrrrr@{}}\\toprule\n" + T("cell type", "tipo celular") + f" & {CYTOKINE_A} & {CYTOKINE_B} & " + T("both", "ambas") + f" & {T('only', 'solo')} {CYTOKINE_A} & {T('only', 'solo')} {CYTOKINE_B} & {CYTOKINE_A} vs {CYTOKINE_B}\\\\\\midrule\n")
    for ct in CELL_TYPES_COMPARE:
        A_, B_, AB_ = cmp[ct]["A"], cmp[ct]["B"], cmp[ct]["AB"]; ok_ = A_["tested"] & B_["tested"]; a_, b_ = (A_["direction"] != "ns") & ok_, (B_["direction"] != "ns") & ok_
        f.write(f"{CT_LABEL[ct]} & {tex_int(a_.sum())} & {tex_int(b_.sum())} & {tex_int((a_ & b_).sum())} & {tex_int((a_ & ~b_).sum())} & {tex_int((~a_ & b_).sum())} & {tex_int((AB_['direction'] != 'ns').sum())}\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n")
print("✅ figures and tex tables written to", HERE)
