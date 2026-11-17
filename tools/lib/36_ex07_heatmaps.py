# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  Exercise 7: heatmaps of raw counts vs normalised values
# ══════════════════════════════════════════════════════════════════════════════════════════════════
HEAT_VIEWS = [("Raw counts", "raw"), ("Normalised: log₂(CPM + 1)", "norm"), ("Both, side by side", "both")]
HEAT_COLS = [("every donor (12 columns per condition)", "donors"), ("the mean of the 12 donors", "mean")]
HEAT_DEFAULT_GROUPS = ["Monocyte markers", "T-cell markers", "CD8 T-cell markers", "NK-cell markers", "IFN-γ / STAT1 response",
                       "IL-4 / STAT6 response", "IL-1β / NF-κB response", "IL-2 / STAT5 response"]
HEAT_DEFAULT_CYTOKINES = ["PBS", "IFN-gamma", "IL-4", "IL-1-beta", "IL-2"]
MAX_HEAT_CONDITIONS, MAX_HEAT_GENES = 25, 70


def zscore_rows(M):
    """z-score every ROW (gene) across its columns (samples): (x − mean) / standard deviation."""
    sd = M.std(axis=1, ddof=1, keepdims=True)
    z = (M - M.mean(axis=1, keepdims=True)) / np.where(sd > 0, sd, 1)
    return np.where(sd > 0, z, 0.0)


def heatmap_data(genes, cell_types, conditions, kind, cols_mode):
    """genes × columns matrix (+ library size per column). Columns are ordered: cell type > condition > donor."""
    cols = [f"{d}|{ct}|{cy}" for ct in cell_types for cy in conditions for d in DONORS]
    idx = adata.obs_names.get_indexer(cols)
    gi = [gene_index(g) for g in genes]
    src = adata.X if kind == "raw" else adata.layers["log2cpm"]
    M = np.asarray(src[np.ix_(idx, gi)], dtype=np.float64).T
    lib = np.log10(adata.obs["library_size"].values[idx].astype(float))
    if cols_mode == "mean":
        M = M.reshape(len(genes), -1, len(DONORS)).mean(axis=2)
        lib = lib.reshape(-1, len(DONORS)).mean(axis=1)
    return M, lib[None, :]


def draw_heatmaps(panels, genes, group_sizes, group_names, cell_types, conditions, cols_mode, lib, title):
    """`panels` = [(title, matrix), …]: same genes and columns, one panel each, with a cell-type bar and a library-size bar on top."""
    n_panel = len(panels)
    n_per_ct = len(conditions) * (len(DONORS) if cols_mode == "donors" else 1)
    n_cols = n_per_ct * len(cell_types)
    height = 0.22 * len(genes) + 2.3
    fig = plt.figure(figsize=(6.4 * n_panel + 1.6, height))
    gs = fig.add_gridspec(3, n_panel + 1, height_ratios=[0.7, 0.7, max(len(genes), 4)], width_ratios=[1] * n_panel + [0.025],
                          hspace=0.06, wspace=0.05, top=1 - 0.55 / height, bottom=0.95 / height + 0.04, left=0.12 if n_panel == 1 else 0.075, right=0.94)
    im = None
    for c, (name, M) in enumerate(panels):
        ax_ct, ax_lib, ax = (fig.add_subplot(gs[r, c]) for r in range(3))
        for k, ct in enumerate(cell_types):
            ax_ct.add_patch(plt.Rectangle((k * n_per_ct - 0.5, 0), n_per_ct, 1, color=CT_COLORS[ct], lw=0))
            ax_ct.text(k * n_per_ct + n_per_ct / 2 - 0.5, 0.5, CT_SHORT[ct], ha="center", va="center", fontsize=9, fontweight="bold", color="black" if ct == "NK_CD56low" else "white")
        ax_ct.set_xlim(-0.5, n_cols - 0.5); ax_ct.set_ylim(0, 1); ax_ct.axis("off"); ax_ct.set_title(name, fontsize=11, pad=6)
        ax_lib.imshow(lib, aspect="auto", cmap="Greys", interpolation="nearest")
        ax_lib.set_xticks([]); ax_lib.set_yticks([0] if c == 0 else [])
        if c == 0:
            ax_lib.set_yticklabels(["library size"], fontsize=7.5)
        for sp in ax_lib.spines.values():
            sp.set_visible(False)
        im = ax.imshow(zscore_rows(M), aspect="auto", cmap=DIVERGING, vmin=-2.5, vmax=2.5, interpolation="nearest")
        ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes if c == 0 else [], fontsize=8, fontstyle="italic")
        if cols_mode == "donors":
            ax.set_xticks([(k + 0.5) * len(DONORS) - 0.5 for k in range(len(cell_types) * len(conditions))])
            ax.set_xticklabels([cyt_label(cy) for _ in cell_types for cy in conditions], rotation=90, fontsize=7)
        else:
            ax.set_xticks(range(n_cols)); ax.set_xticklabels([cyt_label(cy) for _ in cell_types for cy in conditions], rotation=90, fontsize=7.5)
        for b in np.cumsum(group_sizes)[:-1] - 0.5:
            ax.axhline(b, color="white", lw=2)
        for k in range(1, len(cell_types)):
            ax.axvline(k * n_per_ct - 0.5, color="white", lw=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        if c == 0 and group_names:
            starts = np.concatenate([[0], np.cumsum(group_sizes)[:-1]])
            for gname, st, sz in zip(group_names, starts, group_sizes):
                ax.annotate(gname, xy=(0, st + sz / 2 - 0.5), xycoords=("axes fraction", "data"), xytext=(-52, 0), textcoords="offset points",
                            ha="right", va="center", fontsize=7.5, color="#555", annotation_clip=False)
    fig.colorbar(im, cax=fig.add_subplot(gs[2, n_panel]), label="z-score of the gene")
    fig.suptitle(title, y=1 - 0.06 / height, fontsize=10.5, fontweight="bold")
    show_fig(fig)


def ex7_render(view, groups, extra_text, cell_types, cytokines, cols_mode):
    conditions = [c for c in adata.obs["cytokine"].cat.categories if c in set(cytokines)]
    cell_types = [c for c in CELL_TYPES if c in set(cell_types)]
    genes, group_sizes, group_names = [], [], []
    for gname in groups:
        gg = [g for g in GENESETS.get(gname, []) if g not in genes]
        if gg:
            genes += gg; group_sizes.append(len(gg)); group_names.append(gname)
    extra, unknown = parse_genes(extra_text, 40)
    unknown_genes_message(unknown)
    extra = [g for g in extra if g not in genes]
    if extra:
        genes += extra; group_sizes.append(len(extra)); group_names.append("your genes")
    if not genes:
        msg("Select at least one <b>gene group</b> or type a few gene symbols.", "info"); return
    if not conditions or not cell_types:
        msg("Select at least one <b>cell type</b> and one <b>cytokine / condition</b>.", "info"); return
    if len(conditions) > MAX_HEAT_CONDITIONS:
        msg(f"Too many conditions for a readable heatmap: showing the first {MAX_HEAT_CONDITIONS}.", "warn"); conditions = conditions[:MAX_HEAT_CONDITIONS]
    if len(genes) > MAX_HEAT_GENES:
        msg(f"Too many genes for a readable heatmap: showing the first {MAX_HEAT_GENES}.", "warn")
        genes = genes[:MAX_HEAT_GENES]; group_sizes = list(np.minimum(group_sizes, np.maximum(0, MAX_HEAT_GENES - np.concatenate([[0], np.cumsum(group_sizes)[:-1]]))).astype(int))
        group_names = [n for n, s in zip(group_names, group_sizes) if s > 0]; group_sizes = [s for s in group_sizes if s > 0]
    panels, lib = [], None
    for kind, name in (("raw", "BEFORE: raw counts"), ("norm", "AFTER: log₂(CPM + 1)")):
        if view == "both" or view == kind:
            M, lib = heatmap_data(genes, cell_types, conditions, kind, cols_mode)
            panels.append((name, M))
    n_cols = len(cell_types) * len(conditions) * (len(DONORS) if cols_mode == "donors" else 1)
    draw_heatmaps(panels, genes, group_sizes, group_names, cell_types, conditions, cols_mode, lib,
                  f"{len(genes)} genes × {n_cols} columns ({len(cell_types)} cell types × {len(conditions)} conditions" + (" × 12 donors)" if cols_mode == "donors" else ", mean of 12 donors)"))


def exercise_heatmaps():
    cyts = list(adata.obs["cytokine"].cat.categories)
    c = dict(view=widgets.ToggleButtons(options=HEAT_VIEWS, value="both", style={"button_width": "250px"}),
             groups=multi_select(list(GENESETS), pick_default(list(GENESETS), HEAT_DEFAULT_GROUPS), rows=10, width="240px"),
             extra_text=None, cell_types=multi_select([(CT_LABEL[x], x) for x in CELL_TYPES], CELL_TYPES, rows=4, width="170px"),
             cytokines=multi_select(cyts, pick_default(cyts, HEAT_DEFAULT_CYTOKINES), rows=10, width="170px"),
             cols_mode=dropdown(HEAT_COLS, "donors", "290px"))
    ttext, twidget = gene_picker("", "260px", "e.g. TNF, IL1B", with_sets=False, label="Extra genes (optional)")
    c["extra_text"] = ttext
    rows = [[c["view"]],
            [widgets.VBox([field("Gene groups", c["groups"], "Ctrl/⌘ + click: several"), selection_summary(c["groups"], max_items=3)]), widgets.VBox([twidget]),
             widgets.VBox([field("Cell types", c["cell_types"]), field("Columns show", c["cols_mode"])]),
             widgets.VBox([field("Cytokines / conditions", c["cytokines"], "PBS = control"), selection_summary(c["cytokines"], cyt_label)])]]
    ex = Exercise(c, ex7_render, rows, title="What do you want to see?")
    ex.show()
    return ex
