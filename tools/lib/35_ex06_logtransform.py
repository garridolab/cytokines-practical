# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  Exercise 6: logarithms: scale, spread and the zeros
# ══════════════════════════════════════════════════════════════════════════════════════════════════
BASES = [("2  (log₂, the one we use)", "2"), ("e  (natural log)", "e"), ("10  (log₁₀)", "10")]
BASE_VALUE = {"2": 2.0, "e": math.e, "10": 10.0}
PSEUDO_OPTIONS = [0, 0.1, 0.5, 1, 2, 5, 10]
_MV_CACHE = {}


def log_transform(cpm, base, pseudo):
    """log_base(cpm + pseudocount). With pseudocount 0 the zeros become NaN (log of 0 does not exist)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.log(cpm.astype(np.float64) + pseudo) / math.log(BASE_VALUE[base])
    out[~np.isfinite(out)] = np.nan
    return out


def mean_var(M):
    """Mean and variance of every column of M ignoring NaN (columns with fewer than 20 valid values → NaN)."""
    valid = np.isfinite(M).sum(axis=0)
    mean, var = np.nanmean(M, axis=0), np.nanvar(M, axis=0, ddof=1)
    bad = valid < 20
    mean[bad] = np.nan; var[bad] = np.nan
    return mean, var


def top20_share(var):
    v = var[np.isfinite(var)]
    return np.sort(v)[-20:].sum() / v.sum()


def ex6_render(gene, cell_type, base, pseudo):
    g = one_gene(gene, None)
    if g is None:
        msg("Type a valid <b>gene symbol</b> (for example <i>CXCL10</i>, <i>IDO1</i>, <i>ACTB</i>).", "info"); return
    m = adata.obs["cell_type"].values == cell_type
    cpm_all = adata.layers["cpm"][m][:, EXPRESSED]
    key = ("cpm", cell_type)
    if key not in _MV_CACHE:                                                  # the CPM statistics never change: compute them once
        _MV_CACHE[key] = mean_var(cpm_all.astype(np.float64))
    mean_c, var_c = _MV_CACHE[key]
    L = log_transform(cpm_all, base, pseudo)
    mean_l, var_l = mean_var(L)
    j = gene_index(g)
    x = adata.layers["cpm"][m, j].astype(float)
    lx = log_transform(x, base, pseudo)
    n_zero = int((x == 0).sum())
    bname = {"2": "log₂", "e": "ln", "10": "log₁₀"}[base]
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.6))
    col = CT_COLORS[cell_type]
    # (a) the gene, in CPM
    ax = axes[0, 0]
    ax.hist(x, bins=40, color=col, alpha=0.85)
    ax.axvline(np.mean(x), color="k", ls="--", lw=1.5); ax.axvline(np.median(x), color=UP_COLOR, ls=":", lw=2)
    ax.set_xlabel("CPM"); ax.set_ylabel("number of samples")
    sk = stats.skew(x)
    ax.set_title(f"① {g} in CPM: skewed (skewness {sk:.1f})")
    ax.legend(handles=[Line2D([], [], color="k", ls="--", label=f"mean {np.mean(x):,.0f}"), Line2D([], [], color=UP_COLOR, ls=":", lw=2, label=f"median {np.median(x):,.0f}")], loc="upper right")
    # (b) the same gene after the logarithm
    ax = axes[0, 1]
    ok = np.isfinite(lx)
    ax.hist(lx[ok], bins=40, color=col, alpha=0.85)
    ax.set_xlabel(f"{bname}(CPM + {pseudo:g})" if pseudo else f"{bname}(CPM)"); ax.set_ylabel("number of samples")
    sk2 = stats.skew(lx[ok]) if ok.sum() > 2 else float("nan")
    ax.set_title(f"② {g} after the logarithm (skewness {sk2:.1f})")
    ax.text(0.02, 0.95, f"variance of the values:\n{np.var(x, ddof=1):,.0f}  →  {np.var(lx[ok], ddof=1):.2f}", transform=ax.transAxes, ha="left", va="top", fontsize=9.5)
    # (c) all genes: variance vs mean, CPM
    ax = axes[1, 0]
    ok = (mean_c > 0.01) & (var_c > 0)
    ax.scatter(mean_c[ok], var_c[ok], s=3, alpha=0.25, color=BLUE, linewidths=0, rasterized=True)
    xx = np.array([0.01, mean_c[ok].max()]); ax.plot(xx, xx ** 2 * 0.3, color="k", ls="--", lw=1.1, label="variance ∝ mean²")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.legend(loc="upper left")
    ax.set_xlabel("mean CPM of the gene"); ax.set_ylabel("variance across samples")
    ax.set_title(f"③ CPM: variance rises with the mean\n20 genes hold {top20_share(var_c):.0%} of the total variance")
    # (d) all genes: variance vs mean, log
    ax = axes[1, 1]
    ok = np.isfinite(mean_l) & np.isfinite(var_l) & (var_l > 0)
    ax.scatter(mean_l[ok], var_l[ok], s=3, alpha=0.25, color=BLUE, linewidths=0, rasterized=True)
    ax.set_xlabel(f"mean of {bname}(CPM{' + ' + format(pseudo, 'g') if pseudo else ''})"); ax.set_ylabel("variance across samples")
    ax.set_title(f"④ Log: variance not tied to expression\n20 genes hold {top20_share(var_l):.1%} of the total variance")
    fig.suptitle(f"{CT_LABEL[cell_type]} · {m.sum():,} samples · {cpm_all.shape[1]:,} expressed genes", fontsize=10.5, y=1.0)
    fig.tight_layout()
    show_fig(fig)
    notes = []
    if pseudo == 0:
        notes.append(f"With <b>no pseudocount</b>, log(0) does not exist: " + (f"the {n_zero} zero values of {g} ({n_zero / len(x):.0%}) cannot be plotted, and " if n_zero else f"{g} has no zeros here, but ") + "genes with many zeros are dropped from panel ④.")
    else:
        notes.append(f"The pseudocount ({pseudo:g}) is added <i>before</i> the logarithm so that zeros become log({pseudo:g}) instead of −∞. Larger pseudocounts squash the low values even more.")
    if base != "2":
        notes.append(f"Changing the base only <b>rescales</b> the axis (by a constant: log₂(x) = {math.log2(BASE_VALUE[base]):.2f} × {bname}(x)): the shape of the histogram does not change, the variances are multiplied by {(1 / math.log2(BASE_VALUE[base])) ** 2:.2f}).")
    display(HTML('<div class="cx-info">' + "<br>".join(notes) + "</div>"))


def exercise_logtransform():
    c = dict(gene=gene_combo("CXCL10", "170px"), cell_type=dropdown([(CT_LABEL[x], x) for x in CELL_TYPES], "CD14_Mono", "190px"),
             base=dropdown(BASES, "2", "220px"),
             pseudo=widgets.SelectionSlider(options=PSEUDO_OPTIONS, value=1, continuous_update=False, layout=_lay("330px"), style=STYLE))
    rows = [[field("Gene (top panels)", c["gene"], "e.g. CXCL10, IDO1, ACTB"), field("Cell type", c["cell_type"]), field("Base of the logarithm", c["base"]),
             field("Pseudocount added before the log", c["pseudo"], "0 = no pseudocount: log(0) is undefined")]]
    ex = Exercise(c, ex6_render, rows, title="Transform the CPM values")
    ex.show()
    return ex
