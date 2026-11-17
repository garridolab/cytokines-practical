# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  Exercise 5: raw counts vs counts per million (CPM)
# ══════════════════════════════════════════════════════════════════════════════════════════════════
VIEWS5 = [("Real data: same donors, PBS vs a cytokine", "pair"), ("Real data: all samples of a cell type", "all"), ("Toy example (simulated)", "toy")]


def _geo_ratio(a, b):
    """Typical (geometric-mean) ratio b / a over donors, with +1 so that zeros do not break it."""
    return float(np.exp(np.mean(np.log((np.asarray(b, float) + 1) / (np.asarray(a, float) + 1)))))


def _ratio_text(r, what):
    if r >= 1.15:
        return f'<b>{r:.1f}× higher</b>'
    if r <= 1 / 1.15:
        return f'<b>{1 / r:.1f}× lower</b>'
    return '<b>about the same</b>'


def _paired(ax, a, b, labels, title, colors=(PBS_COLOR, ORANGE), log=True):
    for x, y in zip(a, b):
        ax.plot([0, 1], [x, y], color=UP_COLOR if y > x else DOWN_COLOR, lw=1.1, alpha=0.65, zorder=1)
    ax.scatter(np.zeros(len(a)), a, s=26, color=colors[0], zorder=2)
    ax.scatter(np.ones(len(b)), b, s=26, color=colors[1], zorder=2)
    ax.plot([0, 1], [np.exp(np.mean(np.log(np.maximum(a, 0.5)))), np.exp(np.mean(np.log(np.maximum(b, 0.5))))], color="k", lw=3, marker="D", ms=7, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(labels); ax.set_xlim(-0.4, 1.4)
    if log:
        ax.set_yscale("log")
    ax.set_title(title)


def ex5_pair(gene, cell_type, cytokine):
    names = lambda c: [f"{d}|{cell_type}|{c}" for d in DONORS]
    r0, r1 = adata.obs_names.get_indexer(names("PBS")), adata.obs_names.get_indexer(names(cytokine))
    j = gene_index(gene)
    lib0, lib1 = adata.obs["library_size"].values[r0].astype(float), adata.obs["library_size"].values[r1].astype(float)
    raw0, raw1 = adata.X[r0, j].astype(float), adata.X[r1, j].astype(float)
    cpm0, cpm1 = adata.layers["cpm"][r0, j].astype(float), adata.layers["cpm"][r1, j].astype(float)
    cl = cyt_label(cytokine)
    display(HTML(
        f'<div class="cx-info">In <b>{CT_LABEL[cell_type]}</b>, with {cl} the library size is {_ratio_text(_geo_ratio(lib0, lib1), "")} than with PBS. '
        f'The raw counts of <i>{gene}</i> are {_ratio_text(_geo_ratio(raw0, raw1), "")}. '
        f'After CPM, <i>{gene}</i> is {_ratio_text(_geo_ratio(cpm0, cpm1), "")} ({cl} / PBS, typical donor).</div>'))
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9))
    lab = ["PBS", cl]
    _paired(axes[0], lib0, lib1, lab, "① Library size of the sample"); axes[0].set_ylabel("total counts of the sample")
    _paired(axes[1], raw0, raw1, lab, f"② RAW counts of {gene}"); axes[1].set_ylabel("raw counts")
    _paired(axes[2], cpm0, cpm1, lab, f"③ CPM of {gene}  (normalised)"); axes[2].set_ylabel("counts per million (CPM)")
    fig.suptitle(f"{CT_LABEL[cell_type]}: one line per donor (red = higher with {cl}, blue = lower).  Black = typical donor", fontsize=10, y=1.0)
    fig.tight_layout()
    show_fig(fig)


def ex5_all(gene, cell_type):
    m = adata.obs["cell_type"].values == cell_type
    j = gene_index(gene)
    lib = adata.obs["library_size"].values[m].astype(float)
    ctrl = (adata.obs["treatment"].values[m] == "control")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), sharex=True)
    for ax, y, name in zip(axes, (adata.X[m, j].astype(float), adata.layers["cpm"][m, j].astype(float)), ("raw counts", "CPM")):
        ax.scatter(lib[~ctrl], y[~ctrl] + 1, s=8, alpha=0.45, color=CT_COLORS[cell_type], linewidths=0, label="cytokine-treated samples")
        ax.scatter(lib[ctrl], y[ctrl] + 1, s=36, facecolors="none", edgecolors="k", linewidths=1.1, label="PBS controls")
        ax.set_xscale("log"); ax.set_yscale("log")
        ok = y >= 0
        r = np.corrcoef(np.log10(lib), np.log10(y + 1))[0, 1] if np.std(y) > 0 else float("nan")
        if np.isfinite(r):
            slope, icpt = np.polyfit(np.log10(lib), np.log10(y + 1), 1)
            xs = np.array([lib.min(), lib.max()]); ax.plot(xs, 10 ** (icpt + slope * np.log10(xs)), color=UP_COLOR, lw=1.8)
        ax.set_xlabel("library size of the sample (total counts)"); ax.set_ylabel(f"{name} + 1")
        ax.set_title(f"{gene}: {name}\ncorrelation with library size: r = {r:+.2f}")
    axes[0].legend(loc="upper left", markerscale=1.6)
    fig.suptitle(f"{CT_LABEL[cell_type]}: {m.sum():,} samples (every dot is a pseudobulk sample)", fontsize=10, y=1.0)
    fig.tight_layout()
    show_fig(fig)


def ex5_toy(depth, fold):
    names = ["gene A", "gene B\n(interferon gene)", "gene C", "gene D"]
    base = np.array([1000.0, 50.0, 2500.0, 1450.0])                    # sample 1: PBS, deep (5,000 counts in total)
    p2 = base / base.sum(); p2[1] *= fold; p2 /= p2.sum()               # sample 2: gene B truly goes `fold` times up, the others do not change
    n2 = base.sum() * depth
    raw = np.vstack([base, np.round(p2 * n2)])
    cpm = raw / raw.sum(axis=1, keepdims=True) * 1e6
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9))
    x, w = np.arange(4), 0.38
    for ax, M, title, yl in zip(axes, (raw, cpm), ("RAW counts", "COUNTS PER MILLION (CPM)"), ("counts", "CPM")):
        ax.bar(x - w / 2, M[0], w, color=PBS_COLOR, label=f"sample 1: PBS ({raw[0].sum():,.0f} counts in total)")
        ax.bar(x + w / 2, M[1], w, color=ORANGE, label=f"sample 2: with cytokine ({raw[1].sum():,.0f} counts in total)")
        for xi, (a, b) in enumerate(zip(M[0], M[1])):
            ax.text(xi, max(a, b) * 1.03, f"×{b / a:.2f}", ha="center", fontsize=9, fontweight="bold", color=UP_COLOR if b / a > 1.15 else (DOWN_COLOR if b / a < 0.87 else "#333"))
        ax.set_xticks(x); ax.set_xticklabels(names, fontsize=9); ax.set_ylabel(yl); ax.set_title(title)
        ax.set_ylim(0, M.max() * 1.18)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.04), fontsize=9)
    fig.suptitle(f"Only gene B really changes (×{fold:g}). Numbers above the bars: sample 2 / sample 1", fontsize=10, y=1.0)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    show_fig(fig)
    display(HTML(f'<div class="cx-info">Sample 2 was sequenced <b>{depth:.0%}</b> as deeply as sample 1. In raw counts <b>every gene looks '
                 f'{"lower" if depth < 0.95 else "similar"}</b>, even gene B (which truly goes up ×{fold:g}). '
                 f'After dividing by the library size (CPM), genes A, C and D are about unchanged and gene B shows its true increase.</div>'))


def ex5_render(view, gene, cell_type, cytokine, depth, fold):
    if view == "toy":
        return ex5_toy(depth, fold)
    g = one_gene(gene, None)
    if g is None:
        msg("Type a valid <b>gene symbol</b> (for example <i>ACTB</i>, <i>STAT1</i>, <i>IDO1</i>).", "info"); return
    if view == "pair":
        if cytokine == "PBS":
            msg("Choose a cytokine different from PBS.", "info"); return
        ex5_pair(g, cell_type, cytokine)
    else:
        ex5_all(g, cell_type)


def exercise_normalization():
    cyts = [c for c in adata.obs["cytokine"].cat.categories if c != "PBS"]
    c = dict(view=widgets.ToggleButtons(options=VIEWS5, value="pair", style={"button_width": "300px"}, layout=_lay("930px")),
             gene=gene_combo("ACTB", "170px"), cell_type=dropdown([(CT_LABEL[x], x) for x in CELL_TYPES], "CD14_Mono", "190px"),
             cytokine=dropdown([(cyt_label(x), x) for x in cyts], "IFN-gamma", "190px"),
             depth=float_slider(0.1, 1.0, 0.2, 0.05, "280px", ".0%"), fold=float_slider(1, 8, 4, 0.5, "280px", ".1f"))
    f_gene, f_ct, f_cyt = field("Gene", c["gene"], "type a symbol, e.g. ACTB, STAT1, IDO1"), field("Cell type", c["cell_type"]), field("Compare PBS with…", c["cytokine"])
    f_depth = field("Sequencing depth of sample 2 (relative to sample 1)", c["depth"])
    f_fold = field("True fold change of gene B", c["fold"])

    def _toggle(*_):
        real = c["view"].value != "toy"
        for f in (f_gene, f_ct):
            f.layout.display = None if real else "none"
        f_cyt.layout.display = None if c["view"].value == "pair" else "none"
        for f in (f_depth, f_fold):
            f.layout.display = "none" if real else None
    c["view"].observe(_toggle, "value"); _toggle()
    ex = Exercise(c, ex5_render, [[c["view"]], [f_gene, f_ct, f_cyt, f_depth, f_fold]], title="Raw counts or CPM?")
    ex.show()
    return ex
