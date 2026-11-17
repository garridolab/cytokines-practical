# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  Exercise 12: comparing two DE signatures (correlation of log2 fold changes)
# ══════════════════════════════════════════════════════════════════════════════════════════════════
GENE_SUBSETS = [("all genes tested in both comparisons", "all"), ("genes significant in A or B", "either"), ("genes significant in both A and B", "both")]
EXAMPLES12 = [("— choose an example —", ""),
              ("IL-4 vs IL-13 in monocytes: almost the same response", ("CD14_Mono", "IL-4", "CD14_Mono", "IL-13")),
              ("IL-4 vs IL-13 in CD4 T cells: much less similar", ("CD4_T_cell", "IL-4", "CD4_T_cell", "IL-13")),
              ("IFN-β vs IFN-γ in monocytes: two interferons", ("CD14_Mono", "IFN-beta", "CD14_Mono", "IFN-gamma")),
              ("IFN-γ vs IL-4 in monocytes: opposite programmes", ("CD14_Mono", "IFN-gamma", "CD14_Mono", "IL-4")),
              ("TNF-α in monocytes vs in CD4 T cells: same stimulus, two cell types", ("CD14_Mono", "TNF-alpha", "CD4_T_cell", "TNF-alpha")),
              ("IL-6 in monocytes vs in CD4 T cells: same stimulus, two cell types", ("CD14_Mono", "IL-6", "CD4_T_cell", "IL-6")),
              ("TNF-α vs IL-6 in monocytes", ("CD14_Mono", "TNF-alpha", "CD14_Mono", "IL-6")),
              ("IL-2 vs IL-15 in NK cells: the same receptor", ("NK_CD56low", "IL-2", "NK_CD56low", "IL-15"))]
COL_A, COL_B, COL_BOTH = BLUE, ORANGE, PURPLE


def compare_tables(da, db, subset, padj_thr, lfc_thr):
    """Merge two DE tables on the genes tested in both and compute the statistics of the comparison."""
    ok = da["tested"] & db["tested"]
    m = pd.DataFrame({"a": da.loc[ok, "log2fc"], "b": db.loc[ok, "log2fc"], "dir_a": da.loc[ok, "direction"], "dir_b": db.loc[ok, "direction"]})
    sa, sb = m["dir_a"] != "ns", m["dir_b"] != "ns"
    m["cls"] = np.select([sa & sb, sa, sb], ["both", "A", "B"], "none")
    use = {"all": np.ones(len(m), bool), "either": (sa | sb).values, "both": (sa & sb).values}[subset]
    return m, use


def correlations(m, use):
    x, y = m.loc[use, "a"].values.astype(float), m.loc[use, "b"].values.astype(float)
    if len(x) < 5 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan"), float("nan"), len(x)
    return float(np.corrcoef(x, y)[0, 1]), float(stats.spearmanr(x, y)[0]), len(x)


def similarity_ranking(ct, cyt, subset_all=True):
    """Pearson correlation of the log2FC of (ct, cyt) with every other cytokine in the same cell type (cache only)."""
    i, j = DE.cell_types.index(ct), DE.cytokines.index(cyt)
    lfc, nlp = DE.arrays["lfc"][i].astype(np.float32), DE.arrays["nlp_adj"][i].astype(np.float32)
    out = {}
    for k, name in enumerate(DE.cytokines):
        if k == j:
            continue
        ok = ~np.isnan(nlp[j]) & ~np.isnan(nlp[k])
        if ok.sum() > 50:
            out[name] = float(np.corrcoef(lfc[j][ok], lfc[k][ok])[0, 1])
    return pd.Series(out).sort_values(ascending=False)


def ex12_render(cell_a, cyt_a, cell_b, cyt_b, subset, padj_thr, lfc_thr, genes_text, label_top):
    if cyt_a == cyt_b and cell_a == cell_b:
        msg("A and B are the <b>same comparison</b>: the correlation is trivially 1. Change the cytokine or the cell type of B.", "info")
    if not DE.cached and any((c, y) not in DE._memo for c, y in ((cell_a, cyt_a), (cell_b, cyt_b))):
        msg("⏳ The pre-computed results were not found: running DESeq2 (≈ 15 seconds per comparison)…", "info")
    da, db = de_summary(DE.get(cell_a, cyt_a), padj_thr, lfc_thr), de_summary(DE.get(cell_b, cyt_b), padj_thr, lfc_thr)
    m, use = compare_tables(da, db, subset, padj_thr, lfc_thr)
    highlight, unknown = parse_genes(genes_text, 25)
    unknown_genes_message(unknown)
    r, rho, n = correlations(m, use)
    na, nb = f"{cyt_label(cyt_a)} · {CT_SHORT[cell_a]}", f"{cyt_label(cyt_b)} · {CT_SHORT[cell_b]}"
    ranking = similarity_ranking(cell_a, cyt_a) if DE.cached else None
    fig = plt.figure(figsize=(13, 6.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1], wspace=0.35)
    ax = fig.add_subplot(gs[0])
    for cls, col, lab, s_ in (("none", NS_COLOR, "neither", 5), ("A", COL_A, "A only", 11), ("B", COL_B, "B only", 11), ("both", COL_BOTH, "both", 13)):
        mm = m[m["cls"] == cls]
        ax.scatter(mm["a"], mm["b"], s=s_, alpha=0.55 if cls == "none" else 0.8, color=col, linewidths=0, label=f"{disp(lab)}  ({len(mm):,})", rasterized=True, zorder=2)
    ax.axhline(0, color="k", lw=0.6, zorder=1); ax.axvline(0, color="k", lw=0.6, zorder=1)
    lim = float(np.nanpercentile(np.abs(np.r_[m["a"].values, m["b"].values]), 99.8)) * 1.1
    ax.plot([-lim, lim], [-lim, lim], color="#999", lw=0.9, ls="--", zorder=1)
    sel = m[use]
    if len(sel) > 5 and np.std(sel["a"]) > 0:
        slope, icpt = np.polyfit(sel["a"].astype(float), sel["b"].astype(float), 1)
        xs = np.array([-lim, lim]); ax.plot(xs, icpt + slope * xs, color="k", lw=1.8, zorder=4)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    texts = []
    if label_top and len(sel):
        mn = np.minimum(sel["a"].abs(), sel["b"].abs())
        same, opp = sel[np.sign(sel["a"]) == np.sign(sel["b"])], sel[np.sign(sel["a"]) != np.sign(sel["b"])]
        for part, k in ((same, 5), (opp, 3)):
            for g in mn[part.index].nlargest(k).index:
                if g not in highlight:
                    texts.append(ax.text(m.loc[g, "a"], m.loc[g, "b"], g, fontsize=7.5, fontstyle="italic", zorder=5))
    hl = [g for g in highlight if g in m.index]
    if hl:
        ax.scatter(m.loc[hl, "a"], m.loc[hl, "b"], s=70, facecolors="none", edgecolors="k", linewidths=1.6, zorder=6)
        for g in hl:
            texts.append(ax.text(m.loc[g, "a"], m.loc[g, "b"], g, fontsize=9.5, fontweight="bold", fontstyle="italic", zorder=7, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="k", lw=0.6, alpha=0.9)))
    if HAS_ADJUSTTEXT and texts:
        try:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#888", lw=0.5), expand=(1.2, 1.4), max_move=(30, 30))
        except Exception:
            pass
    ax.legend(loc="lower right", markerscale=2.2, fontsize=8.5, title="significant in…", title_fontsize=8.5)
    ax.set_xlabel(f"log₂ fold change of A:  {cyt_label(cyt_a)} vs PBS in {CT_LABEL[cell_a]}"); ax.set_ylabel(f"log₂ fold change of B:  {cyt_label(cyt_b)} vs PBS in {CT_LABEL[cell_b]}")
    ax.set_title(f"Pearson r = {r:.2f}   ·   Spearman ρ = {rho:.2f}   ·   {n:,} genes", fontsize=12)
    ax2 = fig.add_subplot(gs[1])
    if ranking is not None and cell_a == cell_b and cyt_b in ranking.index:
        sel_r = pd.concat([ranking.head(8), ranking.tail(4)])
        if cyt_b not in sel_r.index:
            sel_r = pd.concat([ranking.head(8), ranking.loc[[cyt_b]], ranking.tail(3)])
        cols = [COL_B if k == cyt_b else (UP_COLOR if v > 0 else DOWN_COLOR) for k, v in sel_r.items()]
        ax2.barh([cyt_label(k) for k in sel_r.index], sel_r.values, color=cols)
        ax2.axvline(0, color="k", lw=0.6); ax2.invert_yaxis()
        ax2.set_xlabel("correlation of the log₂ fold changes with A"); ax2.set_title(f"In {CT_LABEL[cell_a]}, who responds like {cyt_label(cyt_a)}?\n(B = {cyt_label(cyt_b)} in orange; rank {int((ranking > ranking[cyt_b]).sum()) + 1} of {len(ranking)})", fontsize=10)
        ax2.tick_params(axis="y", labelsize=8.5)
    elif ranking is not None:
        sel_r = pd.concat([ranking.head(8), ranking.tail(4)])
        ax2.barh([cyt_label(k) for k in sel_r.index], sel_r.values, color=[UP_COLOR if v > 0 else DOWN_COLOR for v in sel_r.values])
        ax2.axvline(0, color="k", lw=0.6); ax2.invert_yaxis(); ax2.tick_params(axis="y", labelsize=8.5)
        ax2.set_xlabel("correlation of the log₂ fold changes with A"); ax2.set_title(f"In {CT_LABEL[cell_a]}, who responds like {cyt_label(cyt_a)}?", fontsize=10)
    else:
        ax2.axis("off")
    fig.tight_layout()
    show_fig(fig)
    # which genes explain the (dis)similarity?
    mn = np.minimum(sel["a"].abs(), sel["b"].abs()) if len(sel) else pd.Series(dtype=float)
    def gene_tab(part, k=8):
        top = mn[part.index].nlargest(k).index
        df = pd.DataFrame({"log₂FC in A": m.loc[top, "a"], "log₂FC in B": m.loc[top, "b"]})
        return html_table(df, index_label="gene")
    if len(sel):
        same_up = sel[(sel["a"] > 0) & (sel["b"] > 0)]; same_dn = sel[(sel["a"] < 0) & (sel["b"] < 0)]; opp = sel[np.sign(sel["a"]) != np.sign(sel["b"])]
        display(HTML(
            f'<div class="cx-info">Among the {n:,} genes used: <b>{len(same_up):,}</b> go up in both, <b>{len(same_dn):,}</b> go down in both, '
            f'<b>{len(opp):,}</b> change in opposite directions. Pearson r = <b>{r:.2f}</b> (Spearman ρ = {rho:.2f}).</div>'
            '<div style="display:flex; flex-wrap:wrap; gap:28px">'
            f'<div><div class="cx-label">Up in both (largest shared change)</div>{gene_tab(same_up)}</div>'
            f'<div><div class="cx-label">Down in both</div>{gene_tab(same_dn)}</div>'
            f'<div><div class="cx-label">Opposite directions</div>{gene_tab(opp, 6)}</div></div>'))
    else:
        msg("No genes pass the selected filter: try <i>all genes tested in both comparisons</i>, or relax the thresholds.", "info")
    if hl:
        df = pd.DataFrame({"log₂FC in A": m.loc[hl, "a"], "adj. p in A": [f"{da.loc[g, 'padj']:.1e}" for g in hl], "log₂FC in B": m.loc[hl, "b"], "adj. p in B": [f"{db.loc[g, 'padj']:.1e}" for g in hl]}, index=hl)
        display(HTML('<div class="cx-label" style="margin-top:8px">Your genes</div>' + html_table(df, index_label="gene")))
    missing = [g for g in highlight if g not in m.index]
    if missing:
        msg("Not analysed in both comparisons (outside the immune-gene list, or too few counts): " + ", ".join(missing), "warn")


def exercise_compare():
    cyts = [c for c in adata.obs["cytokine"].cat.categories if c != "PBS"]
    cts = [(CT_LABEL[x], x) for x in CELL_TYPES]
    cyt_opts = [(cyt_label(x), x) for x in cyts]
    c = dict(cell_a=dropdown(cts, "CD14_Mono", "170px"), cyt_a=dropdown(cyt_opts, "IL-4", "170px"), cell_b=dropdown(cts, "CD14_Mono", "170px"), cyt_b=dropdown(cyt_opts, "IL-13", "170px"),
             subset=dropdown(GENE_SUBSETS, "all", "290px"),
             padj_thr=widgets.SelectionSlider(options=PADJ_OPTIONS, value=PADJ_DEFAULT, continuous_update=False, layout=_lay("260px"), style=STYLE),
             lfc_thr=float_slider(0, 3, LFC_DEFAULT, 0.25, "230px", ".2f"), genes_text=None, label_top=checkbox("label the genes that drive the pattern", True, "300px"))
    gtext, gwidget = gene_picker("", "280px", "e.g. STAT1, CCL17", with_sets=False, label="Highlight genes (optional)")
    c["genes_text"] = gtext
    example = dropdown(EXAMPLES12, "", "480px")
    ex = None

    def _load(change):
        v = change["new"]
        if v:
            ex.set_values(cell_a=v[0], cyt_a=v[1], cell_b=v[2], cyt_b=v[3])
            example.value = ""
    example.observe(_load, "value")
    rows = [[field("Load an example", example, "…or choose A and B yourself below")],
            [widgets.VBox([widgets.HTML('<div class="cx-label" style="color:%s">COMPARISON A</div>' % COL_A), widgets.HBox([field("Cell type", c["cell_a"]), field("Cytokine vs PBS", c["cyt_a"])])]),
             widgets.VBox([widgets.HTML('<div class="cx-label" style="color:%s">COMPARISON B</div>' % COL_B), widgets.HBox([field("Cell type", c["cell_b"]), field("Cytokine vs PBS", c["cyt_b"])])])],
            [field("Genes used for the correlation", c["subset"]), field("“Significant” = adjusted p below…", c["padj_thr"]), field("…and |log₂FC| above", c["lfc_thr"]), gwidget, c["label_top"]]]
    ex = Exercise(c, ex12_render, rows, title="Do two stimuli trigger the same genes?")
    ex.show()
    return ex
