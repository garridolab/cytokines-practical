# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  Exercise 11: volcano plots of the pre-computed DESeq2 results (cytokine vs PBS, one cell type)
# ══════════════════════════════════════════════════════════════════════════════════════════════════
PADJ_OPTIONS = [0.001, 0.01, 0.05, 0.1, 0.2]
try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except Exception:
    HAS_ADJUSTTEXT = False


def de_summary(res, padj_thr, lfc_thr):
    """Add the up / down / ns call and the -log10 adjusted p-value to a DE table."""
    d = res.copy()
    d["direction"] = classify_de(d, padj_thr, lfc_thr)
    d["nlp"] = d["neglog10_padj"]
    d["score"] = d["nlp"] * d["log2fc"].abs()                                    # evidence × effect size: used to pick the genes to label
    return d


def mean_expression(gene_list, cell_type, cytokine):
    """Mean log2(CPM + 1) of each gene in the 12 PBS and the 12 cytokine samples of one cell type."""
    r0 = adata.obs_names.get_indexer([f"{d}|{cell_type}|PBS" for d in DONORS])
    r1 = adata.obs_names.get_indexer([f"{d}|{cell_type}|{cytokine}" for d in DONORS])
    gi = [gene_index(g) for g in gene_list]
    L = adata.layers["log2cpm"]
    return L[np.ix_(r0, gi)].mean(axis=0), L[np.ix_(r1, gi)].mean(axis=0)


def axis_caps(x, y):
    """A few genes have extreme values (e.g. log2FC = 12 when a gene is off in PBS): they are drawn on the edge instead of squeezing the plot."""
    xq, yq = np.quantile(np.abs(x), 0.999), np.quantile(y, 0.999)
    return min(np.abs(x).max(), max(2.0 * xq, 3.0)) * 1.04, min(y.max(), max(1.8 * yq, 10.0)) * 1.04


def draw_volcano(ax, d, padj_thr, lfc_thr, n_labels, highlight=(), title=""):
    t = d[d["tested"]].copy()
    xcap, ycap = axis_caps(t["log2fc"].values, t["nlp"].values)
    t["x"], t["y"] = t["log2fc"].clip(-xcap, xcap), t["nlp"].clip(upper=ycap)
    beyond = (t["x"] != t["log2fc"]) | (t["y"] != t["nlp"])
    for k, col, size, alpha in (("ns", NS_COLOR, 6, 0.55), ("down", DOWN_COLOR, 12, 0.8), ("up", UP_COLOR, 12, 0.8)):
        sub = t[(t["direction"] == k) & ~beyond]
        ax.scatter(sub["x"], sub["y"], s=size, color=col, alpha=alpha, linewidths=0, rasterized=True, zorder=2)
        sub = t[(t["direction"] == k) & beyond]
        ax.scatter(sub["x"], sub["y"], s=40, color=col, marker="^", linewidths=0.5, edgecolors="k", zorder=3)
    for x in (-lfc_thr, lfc_thr):
        ax.axvline(x, color="#999", lw=0.8, ls="--", zorder=1)
    ax.axhline(-np.log10(padj_thr), color="#999", lw=0.8, ls="--", zorder=1)
    ax.set_xlim(-xcap * 1.04, xcap * 1.04); ax.set_ylim(-ycap * 0.03, ycap * 1.1)
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    sig = t[t["direction"] != "ns"]
    hl = [g for g in highlight if g in t.index]
    texts = []
    if n_labels > 0 and len(sig):                                               # the strongest DEGs, half up and half down
        for k in ("up", "down"):
            top = sig[sig["direction"] == k].nlargest(max(1, n_labels // 2), "score")
            for g, r in top.iterrows():
                if g not in hl and not (r["nlp"] > y0 + 0.85 * (y1 - y0) and (r["log2fc"] < x0 + 0.25 * (x1 - x0) or r["log2fc"] > x1 - 0.25 * (x1 - x0))):
                    texts.append(ax.text(r["x"], r["y"], g, fontsize=7.5, fontstyle="italic", zorder=5))
    if hl:                                                                       # the genes chosen by the student
        r = t.loc[hl]
        ax.scatter(r["x"], r["y"], s=70, facecolors="none", edgecolors="k", linewidths=1.6, zorder=6)
        for g, rr in r.iterrows():
            texts.append(ax.text(rr["x"], rr["y"], g, fontsize=9.5, fontweight="bold", fontstyle="italic", zorder=7,
                                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="k", lw=0.6, alpha=0.9)))
    if HAS_ADJUSTTEXT and texts:
        try:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#888", lw=0.5), expand=(1.2, 1.4), max_move=(30, 30))
        except Exception:
            pass
    n_up, n_down = int((t["direction"] == "up").sum()), int((t["direction"] == "down").sum())
    ax.annotate(f"↑ {n_up:,} up", xy=(0.98, 0.98), xycoords="axes fraction", ha="right", va="top", color=UP_COLOR, fontsize=11, fontweight="bold")
    ax.annotate(f"↓ {n_down:,} down", xy=(0.02, 0.98), xycoords="axes fraction", ha="left", va="top", color=DOWN_COLOR, fontsize=11, fontweight="bold")
    foot = f"{len(t):,} genes tested · {len(t) - n_up - n_down:,} not significant" + (f"   ·   ▲ {int(beyond.sum())} beyond the axes (see the tables)" if beyond.any() else "")
    ax.text(0.5, 1.005, foot, transform=ax.transAxes, ha="center", va="bottom", color="#777", fontsize=8.5)
    ax.set_xlabel("log₂ fold change  (cytokine / PBS)"); ax.set_ylabel("−log₁₀ adjusted p-value")
    ax.set_title(title, pad=16)
    return n_up, n_down


def _status(r, padj_thr, lfc_thr):
    if not r["tested"]:
        return "not analysed (outside the immune-gene list, or too few counts)"
    return {"up": '<span class="cx-up">↑ up-regulated</span>', "down": '<span class="cx-down">↓ down-regulated</span>', "ns": "not significant"}[r["direction"]]


def gene_stats_table(d, genes, cell_type, cytokine, padj_thr, lfc_thr):
    m0, m1 = mean_expression(genes, cell_type, cytokine)
    rows = []
    for g, a, b in zip(genes, m0, m1):
        r = d.loc[g]
        rows.append({"gene": g, "log₂ fold change": r["log2fc"] if r["tested"] else np.nan, "fold change": (2 ** r["log2fc"]) if r["tested"] else np.nan,
                     "adjusted p-value": f"{r['padj']:.1e}" if r["tested"] else "—", "mean in PBS": a, f"mean with {cyt_label(cytokine)}": b, "verdict": _status(r, padj_thr, lfc_thr)})
    df = pd.DataFrame(rows).set_index("gene")
    return html_table(df, float_fmt="{:,.2f}", left_cols=["verdict"], index_label=disp("gene"))


def top_table(d, direction, n=10):
    t = d[(d["direction"] == direction)].nlargest(n, "score")
    df = pd.DataFrame({"log₂ fold change": t["log2fc"], "fold change": 2.0 ** t["log2fc"], "adjusted p-value": t["padj"].map(lambda x: f"{x:.1e}")})
    return html_table(df, float_fmt="{:,.2f}", index_label=disp("gene"))


def ex11_render(cell_type, cytokine, padj_thr, lfc_thr, n_labels, genes_text):
    if cytokine not in DE.cytokines and DE.cached:
        msg("Choose a cytokine from the list.", "info"); return
    highlight, unknown = parse_genes(genes_text, 25)
    unknown_genes_message(unknown)
    if not DE.cached and (cell_type, cytokine) not in DE._memo:
        msg("⏳ The pre-computed results were not found: running DESeq2 for this comparison (≈ 15 seconds)…", "info")
    res = DE.get(cell_type, cytokine)
    d = de_summary(res, padj_thr, lfc_thr)
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    draw_volcano(ax, d, padj_thr, lfc_thr, n_labels, highlight, f"{cyt_label(cytokine)} vs PBS · {CT_LABEL[cell_type]}")
    fig.tight_layout()
    show_fig(fig)
    if highlight:
        display(HTML(f'<div class="cx-label">Your genes</div>' + gene_stats_table(d, highlight, cell_type, cytokine, padj_thr, lfc_thr)))
    display(HTML('<div style="display:flex; flex-wrap:wrap; gap:28px"><div><div class="cx-label">Top 10 up-regulated genes <span class="cx-hint">(strongest evidence × effect)</span></div>'
                 + top_table(d, "up") + '</div><div><div class="cx-label">Top 10 down-regulated genes</div>' + top_table(d, "down") + '</div></div>'))


def exercise_volcano():
    cyts = [c for c in adata.obs["cytokine"].cat.categories if c != "PBS"]
    c = dict(cell_type=dropdown([(CT_LABEL[x], x) for x in CELL_TYPES], "CD14_Mono", "190px"), cytokine=dropdown([(cyt_label(x), x) for x in cyts], "IFN-gamma", "190px"),
             padj_thr=widgets.SelectionSlider(options=PADJ_OPTIONS, value=PADJ_DEFAULT, continuous_update=False, layout=_lay("300px"), style=STYLE),
             lfc_thr=float_slider(0, 3, LFC_DEFAULT, 0.25, "270px", ".2f"), n_labels=int_slider(0, 30, 12, 2, "230px"), genes_text=None)
    gtext, gwidget = gene_picker("CXCL9, CXCL10, IDO1, STAT1, IRF1", "300px", "e.g. CXCL9, IDO1", label="Search and highlight genes")
    c["genes_text"] = gtext
    rows = [[field("Cell type", c["cell_type"]), field("Cytokine (vs PBS)", c["cytokine"]), field("Significance: adjusted p-value below…", c["padj_thr"]),
             field("Effect size: |log₂ fold change| above…", c["lfc_thr"], "0.5 ≈ 1.4× · 1 = 2× · 2 = 4×"), field("Gene labels", c["n_labels"], "strongest genes")], [gwidget]]
    ex = Exercise(c, ex11_render, rows, title="Which genes change with a cytokine?")
    ex.show()
    return ex
