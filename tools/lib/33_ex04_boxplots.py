# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  Exercise 4: distributions of raw counts (boxplots of pseudobulk samples)
# ══════════════════════════════════════════════════════════════════════════════════════════════════
SPLITS = [("cell type", "cell_type"), ("control (PBS) vs treated", "treatment"), (disp("donor"), "donor"), ("cytokine (the ones selected below)", "cytokine")]
SPLIT_NAME = {"cell_type": "cell type", "treatment": disp("condition"), "donor": disp("donor"), "cytokine": disp("cytokine")}


def group_colors(split_by, groups):
    if split_by == "cell_type":
        return {g: CT_COLORS[g] for g in groups}
    if split_by == "treatment":
        return {"control": "#8a8984", "treated": ORANGE}
    pal = DONOR_HUES if split_by == "donor" else CYT_HUES
    return {g: pal[i % len(pal)] for i, g in enumerate(groups)}


def group_label(split_by, g):
    return CT_LABEL.get(g, g) if split_by == "cell_type" else (cyt_label(g) if split_by == "cytokine" else disp(g))


def draw_boxes(ax, data, genes, groups, colors, log_y=True, points=True, label_fn=lambda g: g):
    """Boxes of `data[gene][group]` (arrays), grouped by gene on the x axis and coloured by group."""
    G = len(groups)
    w = 0.8 / G
    rng = np.random.default_rng(0)
    floor = 0.5
    for i, gene in enumerate(genes):
        for j, g in enumerate(groups):
            v = np.asarray(data[gene][g], float)
            if len(v) == 0:
                continue
            v = np.maximum(v, floor) if log_y else v
            pos = i + (j - (G - 1) / 2) * w
            bp = ax.boxplot([v], positions=[pos], widths=w * 0.85, patch_artist=True, showfliers=False, medianprops=dict(color="k", lw=1.6),
                            whiskerprops=dict(color="#666", lw=0.9), capprops=dict(color="#666", lw=0.9))
            bp["boxes"][0].set(facecolor=colors[g], alpha=0.55, edgecolor="#555", lw=0.9)
            if points:
                vv = v if len(v) <= 120 else rng.choice(v, 120, replace=False)
                ax.scatter(pos + rng.uniform(-w * 0.3, w * 0.3, len(vv)), vv, s=5, color=colors[g], alpha=0.55, linewidths=0, zorder=3)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, fontsize=9.5, fontstyle="italic")
    ax.set_xlim(-0.6, len(genes) - 0.4)
    if log_y:
        ax.set_yscale("log")


def ex4_render(genes_text, cell_types, split_by, all_cytokines, cytokines, log_y, show_points):
    genes, unknown = parse_genes(genes_text, 12)
    unknown_genes_message(unknown)
    cell_types, cytokines = list(cell_types), list(cytokines)
    if not genes:
        msg("Type at least one valid <b>gene symbol</b> (for example <i>ACTB, LYZ, CD3E</i>).", "info"); return
    if not cell_types:
        msg("Select at least one <b>cell type</b>.", "info"); return
    if split_by == "cytokine" and not cytokines:
        msg("To split by cytokine, select at least one cytokine in the list (hold Ctrl/⌘ to select several).", "info"); return
    m = sample_mask(cell_types, None if (all_cytokines and split_by != "cytokine") else (cytokines or None))
    if m.sum() == 0:
        msg("No samples match this selection.", "warn"); return
    obs = adata.obs[m]
    if split_by == "cytokine":
        groups = [c for c in obs["cytokine"].cat.categories if c in cytokines][:10]
    else:
        groups = [g for g in obs[split_by].cat.categories if g in set(obs[split_by])]
    keep = obs[split_by].isin(groups).values
    rows = np.flatnonzero(m)[keep]
    labels = obs[split_by].astype(str).values[keep]
    colors = group_colors(split_by, groups)
    cols = [gene_index(g) for g in genes]
    block = adata.X[np.ix_(rows, cols)]
    data = {gn: {g: block[labels == g, k] for g in groups} for k, gn in enumerate(genes)}
    lib = {g: adata.obs["library_size"].values[rows][labels == g].astype(float) for g in groups}
    G, nG = len(groups), len(genes)
    fig, axes = plt.subplots(1, 2, figsize=(max(8.5, 1.15 * nG * max(G, 2) ** 0.5 + 4.2), 4.6), gridspec_kw={"width_ratios": [max(nG * (0.5 + 0.3 * G), 3), max(G, 2)]})
    draw_boxes(axes[0], data, genes, groups, colors, log_y, show_points)
    axes[0].set_ylabel("raw counts of the gene in a sample" + ("  (log scale)" if log_y else ""))
    axes[0].set_title("Each box = the pseudobulk samples of one group")
    ax = axes[1]
    for j, g in enumerate(groups):
        v = lib[g]
        bp = ax.boxplot([v], positions=[j], widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color="k", lw=1.6))
        bp["boxes"][0].set(facecolor=colors[g], alpha=0.55, edgecolor="#555")
        if show_points:
            rng = np.random.default_rng(j)
            vv = v if len(v) <= 120 else rng.choice(v, 120, replace=False)
            ax.scatter(j + rng.uniform(-0.2, 0.2, len(vv)), vv, s=5, color=colors[g], alpha=0.55, linewidths=0, zorder=3)
    ax.set_xticks(range(G)); ax.set_xticklabels([group_label(split_by, g) for g in groups], rotation=35 if G > 2 else 0, ha="right" if G > 2 else "center", fontsize=8.5)
    ax.set_yscale("log" if log_y else "linear"); ax.set_ylabel("library size (total counts of the sample)"); ax.set_title("Total counts per sample")
    handles = [Line2D([], [], marker="s", ls="", color=colors[g], alpha=0.7, label=f"{group_label(split_by, g)}  (n = {int((labels == g).sum())})") for g in groups]
    fig.legend(handles=handles, loc="lower center", ncol=min(G, 5), bbox_to_anchor=(0.5, -0.06), fontsize=8.5, title=f"colour = {SPLIT_NAME[split_by]}", title_fontsize=9)
    fig.tight_layout()
    show_fig(fig)
    # the numbers behind the boxes
    tab = pd.DataFrame({group_label(split_by, g): [np.median(data[gn][g]) for gn in genes] + [np.median(lib[g])] for g in groups}, index=genes + ["library size"])
    tab = tab.map(lambda x: f"{x:,.0f}") if hasattr(tab, "map") else tab.applymap(lambda x: f"{x:,.0f}")
    display(HTML('<div class="cx-label">Median value of each box (raw counts)</div>' + html_table(tab, index_label=disp("gene"))))


def exercise_boxplots():
    gtext, gwidget = gene_picker("ACTB, GAPDH, LYZ, CD3E, NKG7, STAT1", placeholder="ACTB, LYZ, ...")
    cyts = list(adata.obs["cytokine"].cat.categories)
    c = dict(genes_text=gtext, cell_types=multi_select([(CT_LABEL[x], x) for x in CELL_TYPES], CELL_TYPES, rows=4, width="170px"),
             split_by=dropdown(SPLITS, "cell_type", "270px"), all_cytokines=checkbox("use all 88 conditions", True, "190px"),
             cytokines=multi_select(cyts, pick_default(cyts, ["PBS", "IFN-gamma", "IL-4", "IL-1-beta"]), rows=6, width="170px"),
             log_y=checkbox("log scale (y axis)", True, "170px"), show_points=checkbox("show the samples as dots", True, "200px"))

    def _sync(*_):
        c["cytokines"].disabled = c["all_cytokines"].value and c["split_by"].value != "cytokine"
    c["all_cytokines"].observe(_sync, "value"); c["split_by"].observe(_sync, "value"); _sync()
    rows = [[widgets.VBox([gwidget]), field("Cell types", c["cell_types"]), field("Group the boxes by", c["split_by"]),
             widgets.VBox([field("Cytokines", c["cytokines"]), selection_summary(c["cytokines"], cyt_label), c["all_cytokines"]])], [c["log_y"], c["show_points"]]]
    ex = Exercise(c, ex4_render, rows, title="Which genes, which samples?")
    ex.show()
    return ex
