# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  Exercises 9 and 10: PCA scatter plots (all samples, and one cell type at a time)
# ══════════════════════════════════════════════════════════════════════════════════════════════════
def pca_scatter(ax, P, color_by, pcs=(0, 1), cytokines=(), gene=None, s=14):
    """Scatter of two PCs. color_by: cell_type | donor | treatment | cytokine | n_cells | library_size | gene."""
    sc, (a, b), obs = P["scores"], pcs, P["obs"]
    kw = dict(linewidths=0.3, edgecolors="white")
    legend = True
    if color_by == "cell_type":
        for ct in CELL_TYPES:
            m = (obs["cell_type"] == ct).values
            if m.any():
                ax.scatter(sc[m, a], sc[m, b], s=s, alpha=0.8, color=CT_COLORS[ct], marker=CT_MARKERS[ct], label=CT_LABEL[ct], **kw)
    elif color_by == "donor":
        for i, dn in enumerate(DONORS):
            m = (obs["donor"] == dn).values
            ax.scatter(sc[m, a], sc[m, b], s=s, alpha=0.85, color=DONOR_HUES[i % 6], marker="o" if i < 6 else "^", label=dn, **kw)
    elif color_by == "treatment":
        m = (obs["treatment"] == "treated").values
        ax.scatter(sc[m, a], sc[m, b], s=s * 0.8, color=NS_COLOR, alpha=0.6, linewidths=0, label=f"treated with a cytokine ({m.sum():,})", zorder=1)
        ax.scatter(sc[~m, a], sc[~m, b], s=s * 2.6, color="black", label=f"control, PBS ({(~m).sum():,})", zorder=2, **kw)
    elif color_by == "cytokine":
        order = {c: i for i, c in enumerate(adata.obs["cytokine"].cat.categories)}
        picked = sorted(cytokines, key=lambda c: order.get(c, 999))[:len(CYT_HUES)]          # same colour for the same cytokine, whatever the click order
        m = ~obs["cytokine"].isin(["PBS"] + picked).values
        ax.scatter(sc[m, a], sc[m, b], s=s * 0.6, color=NS_COLOR, alpha=0.5, linewidths=0, label="all other cytokines", zorder=1)
        for i, cy in enumerate([c for c in picked if c != "PBS"]):
            m = (obs["cytokine"] == cy).values
            ax.scatter(sc[m, a], sc[m, b], s=s * 1.5, color=CYT_HUES[i], marker=CYT_MARKERS[i], label=cyt_label(cy), zorder=2, **kw)
        m = (obs["cytokine"] == "PBS").values
        ax.scatter(sc[m, a], sc[m, b], s=s * 2.6, color="black", label="PBS (control)", zorder=3, **kw)
    else:                                                                       # a number: n_cells, library_size or a gene
        legend = False
        if color_by == "gene":
            values, label = P_gene_values(P, gene), f"{gene}\nlog₂(CPM + 1)"
        elif color_by == "library_size":
            values, label = np.log10(obs["library_size"].values.astype(float)), "log₁₀ library size"
        else:
            values, label = np.log10(obs["psbulk_n_cells"].values.astype(float)), "log₁₀ number of cells"
        order = np.argsort(values)                                              # highest values on top
        pts = ax.scatter(sc[order, a], sc[order, b], c=values[order], cmap=SEQ_POINTS, s=s * 1.5, **kw)
        plt.colorbar(pts, ax=ax, shrink=0.8, pad=0.01, label=label)
    if legend:
        n_leg = len(ax.get_legend_handles_labels()[0])
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=8.5, markerscale=1.5, handletextpad=0.3, columnspacing=1.2, ncol=min(n_leg, 6))
    ax.set_xlabel(f"PC{a + 1}  ({P['var'][a]:.1f} % of the variance)"); ax.set_ylabel(f"PC{b + 1}  ({P['var'][b]:.1f} %)")
    ax.axhline(0, color="#ddd", lw=0.6, zorder=0); ax.axvline(0, color="#ddd", lw=0.6, zorder=0)


def P_gene_values(P, gene):
    return adata.layers["log2cpm"][P["mask"], gene_index(gene)].astype(float)


def plot_scree(ax, P, pcs, n=10):
    var = P["var"][:n]
    ax.bar(range(1, len(var) + 1), var, color=[ORANGE if i in pcs else "#9bbbe6" for i in range(len(var))])
    for i, v in enumerate(var[:6]):
        ax.text(i + 1, v + 0.8, f"{v:.0f}%" if v >= 10 else f"{v:.1f}%", ha="center", fontsize=7.5)
    ax.set_xticks(range(1, len(var) + 1)); ax.set_xlabel("principal component"); ax.set_ylabel("variance explained (%)"); ax.set_title("How much variance does each PC hold?")


def plot_r2(ax, P, pcs, factors=None):
    """What does each of the chosen PCs capture? R² with every sample variable."""
    R2 = P["r2"]
    cols = [c for c in (factors or R2.columns) if c in R2.columns]
    y = np.arange(len(cols))
    w = 0.38
    for k, (pc, col) in enumerate(zip(pcs, (ORANGE, BLUE))):
        vals = [R2.loc[f"PC{pc + 1}", c] for c in cols]
        ax.barh(y + (k - 0.5) * w, vals, w, color=col, label=f"PC{pc + 1}")
        for yy, v in zip(y + (k - 0.5) * w, vals):
            if v > 0.03:
                ax.text(v + 0.01, yy, f"{v:.2f}", va="center", fontsize=7.5)
    ax.set_yticks(y); ax.set_yticklabels([disp(c) for c in cols], fontsize=8.5); ax.invert_yaxis(); ax.set_xlim(0, 1.12)
    ax.set_xlabel("R²  (fraction of the PC explained by the variable)"); ax.set_title("What does each PC capture?"); ax.legend(loc="lower right")


COLOR_CHOICES = [(disp("Donor"), "donor"), ("Cell type", "cell_type"), ("Control vs treated", "treatment"), ("Cytokine (highlight a few)", "cytokine"),
                 ("Number of cells in the sample", "n_cells"), ("Library size", "library_size"), ("Expression of a gene", "gene")]
PC_CHOICES = [(f"PC{i}", i - 1) for i in range(1, 11)]


def _pca_common_msgs(color_by, gene, cytokines):
    if color_by == "gene":
        g = one_gene(gene, None)
        if g is None:
            msg("Type a valid <b>gene symbol</b> to colour by its expression (e.g. <i>LYZ</i>, <i>CD3E</i>, <i>UTY</i>).", "info")
        return g
    if color_by == "cytokine" and len(cytokines) > len(CYT_HUES):
        msg(f"Only the first {len(CYT_HUES)} cytokines are coloured (more colours cannot be told apart).", "warn")
    return None


def ex9_render(color_by, pc_x, pc_y, cytokines, gene):
    g = _pca_common_msgs(color_by, gene, cytokines)
    if color_by == "gene" and g is None:
        return
    if pc_x == pc_y:
        msg("Choose two <b>different</b> PCs for the two axes.", "info"); return
    P = PCA_ALL
    fig = plt.figure(figsize=(13, 5.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1], height_ratios=[1, 1.15], hspace=0.55, wspace=0.28)
    ax = fig.add_subplot(gs[:, 0])
    pca_scatter(ax, P, color_by, (pc_x, pc_y), cytokines, g)
    ax.set_title(f"All {adata.n_obs:,} samples · PCA of the {N_HVG:,} most variable genes")
    plot_scree(fig.add_subplot(gs[0, 1]), P, (pc_x, pc_y))
    plot_r2(fig.add_subplot(gs[1, 1]), P, (pc_x, pc_y), ["cell type", "donor", "cytokine", "control vs treated", "library size", "number of cells"])
    show_fig(fig)


def exercise_pca_all():
    cyts = [c for c in adata.obs["cytokine"].cat.categories if c != "PBS"]
    c = dict(color_by=dropdown(COLOR_CHOICES, "donor", "280px"), pc_x=dropdown(PC_CHOICES, 0, "100px"), pc_y=dropdown(PC_CHOICES, 1, "100px"),
             cytokines=multi_select(cyts, pick_default(cyts, ["IFN-gamma", "IL-4", "TNF-alpha", "IL-10", "IL-6", "GM-CSF"]), rows=7, width="170px"),
             gene=gene_combo("LYZ", "160px"))
    f_cyt = widgets.VBox([field("Cytokines to highlight", c["cytokines"], "up to 8; PBS is always black"), selection_summary(c["cytokines"], cyt_label, 4)])
    f_gene = field("Gene", c["gene"], "type a symbol")

    def _toggle(*_):
        f_cyt.layout.display = None if c["color_by"].value == "cytokine" else "none"
        f_gene.layout.display = None if c["color_by"].value == "gene" else "none"
    c["color_by"].observe(_toggle, "value"); _toggle()
    ex = Exercise(c, ex9_render, [[field("Colour the samples by…", c["color_by"]), field("x axis", c["pc_x"]), field("y axis", c["pc_y"]), f_cyt, f_gene]], title="PCA of the complete dataset")
    ex.show()
    return ex


# ── Exercise 10: one cell type at a time ──────────────────────────────────────────────────────────
COLOR_CHOICES_CT = [("Cytokine (highlight a few)", "cytokine"), (disp("Donor"), "donor"), ("Control vs treated", "treatment"), ("Library size", "library_size"), ("Expression of a gene", "gene")]


def plot_loadings(ax, P, pc, n=8):
    """The genes with the largest weights (loadings) of one PC."""
    s = pd.Series(P["loadings"][:, pc], index=P["hvg_names"]).sort_values()
    top = pd.concat([s.head(n), s.tail(n)])
    ax.barh(top.index, top.values, color=[DOWN_COLOR if v < 0 else UP_COLOR for v in top.values])
    ax.axvline(0, color="k", lw=0.6); ax.tick_params(axis="y", labelsize=8)
    for lab in ax.get_yticklabels():
        lab.set_fontstyle("italic")
    ax.set_xlabel("loading (weight of the gene)"); ax.set_title(f"Genes behind PC{pc + 1}")


def ex10_render(cell_type, color_by, pc_x, pc_y, cytokines, gene, show_loadings):
    g = _pca_common_msgs(color_by, gene, cytokines)
    if color_by == "gene" and g is None:
        return
    if pc_x == pc_y:
        msg("Choose two <b>different</b> PCs for the two axes.", "info"); return
    P = PCA_CT[cell_type]
    if show_loadings:
        fig = plt.figure(figsize=(13.5, 5.4))
        gs = fig.add_gridspec(2, 2, width_ratios=[1.6, 1], hspace=0.55, wspace=0.55)
        ax = fig.add_subplot(gs[:, 0])
        plot_loadings(fig.add_subplot(gs[0, 1]), P, pc_x)
        plot_loadings(fig.add_subplot(gs[1, 1]), P, pc_y)
    else:
        fig, ax = plt.subplots(figsize=(8.6, 5.2))
    pca_scatter(ax, P, color_by, (pc_x, pc_y), cytokines, g, s=22)
    ax.set_title(f"{CT_LABEL[cell_type]} only ({P['obs'].shape[0]:,} samples): PCA of its {N_HVG:,} most variable genes")
    show_fig(fig)
    R2 = P["r2"]
    rows = {f"PC{pc + 1}": {"variance": f"{P['var'][pc]:.1f} %", "R² donor": f"{R2.loc[f'PC{pc + 1}', 'donor']:.2f}", "R² cytokine": f"{R2.loc[f'PC{pc + 1}', 'cytokine']:.2f}",
                              "R² library size": f"{R2.loc[f'PC{pc + 1}', 'library size']:.2f}"} for pc in (pc_x, pc_y)}
    display(HTML('<div class="cx-label">What do these two PCs capture in this cell type?</div>' + html_table(pd.DataFrame(rows).T, index_label="PC")))


def exercise_pca_celltype():
    cyts = [c for c in adata.obs["cytokine"].cat.categories if c != "PBS"]
    c = dict(cell_type=dropdown([(CT_LABEL[x], x) for x in CELL_TYPES], "CD14_Mono", "190px"), color_by=dropdown(COLOR_CHOICES_CT, "cytokine", "250px"),
             pc_x=dropdown(PC_CHOICES, 0, "100px"), pc_y=dropdown(PC_CHOICES, 3, "100px"),
             cytokines=multi_select(cyts, pick_default(cyts, ["IFN-beta", "IFN-gamma", "IL-4", "IL-13", "GM-CSF", "IL-10"]), rows=7, width="170px"),
             gene=gene_combo("UTY", "160px"), show_loadings=checkbox("show the genes behind each axis", True, "250px"))
    f_cyt = widgets.VBox([field("Cytokines to highlight", c["cytokines"], "up to 8; PBS is always black"), selection_summary(c["cytokines"], cyt_label, 4)])
    f_gene = field("Gene", c["gene"], "type a symbol")

    def _toggle(*_):
        f_cyt.layout.display = None if c["color_by"].value == "cytokine" else "none"
        f_gene.layout.display = None if c["color_by"].value == "gene" else "none"
    c["color_by"].observe(_toggle, "value"); _toggle()
    ex = Exercise(c, ex10_render, [[field("Cell type", c["cell_type"]), field("Colour the samples by…", c["color_by"]), field("x axis", c["pc_x"]), field("y axis", c["pc_y"]), f_cyt, f_gene],
                                   [c["show_loadings"]]], title="PCA within one cell type")
    ex.show()
    return ex
