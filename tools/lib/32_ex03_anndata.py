# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  Exercise 3: exploring the AnnData object (X, obs, var) and subsetting it
# ══════════════════════════════════════════════════════════════════════════════════════════════════
class SamplePicker:
    """Pick samples from the 4,224 sample IDs: three dropdowns narrow the list, and the selection survives the filtering."""
    ALL = "(all)"

    def __init__(self, default, rows=9, width="400px", max_options=400):
        obs = adata.obs
        order = np.lexsort((obs["cytokine"].cat.codes.values, obs["cell_type"].cat.codes.values, obs["donor"].cat.codes.values))
        self.ids = list(obs.index[order])                                       # by donor, then cell type, then cytokine
        self.rank = {s: i for i, s in enumerate(self.ids)}
        self.pos = adata.obs_names.get_indexer(self.ids)
        self.max_options = max_options
        self.chosen = [s for s in default if s in self.rank]
        self.f_donor = dropdown([self.ALL] + DONORS, self.ALL, "130px")
        self.f_ct = dropdown([self.ALL] + CELL_TYPES, self.ALL, "150px")
        self.f_cyt = dropdown([self.ALL] + list(obs["cytokine"].cat.categories), self.ALL, "150px")
        self.select = multi_select([], [], rows=rows, width=width)
        self.value_widget = widgets.Text(value="\n".join(self.chosen))             # (never displayed) the selection, as ONE value that only
        self._busy = False                                                          # changes when the student selects, not when the list is filtered
        self.status = widgets.HTML()
        for f in (self.f_donor, self.f_ct, self.f_cyt):
            f.observe(self._refresh, "value")
        self.select.observe(self._on_select, "value")
        self._refresh()
        self.widget = widgets.VBox([
            widgets.HBox([field("Donor", self.f_donor), field("Cell type", self.f_ct), field("Cytokine", self.f_cyt)]),
            self.select, self.status])

    def _refresh(self, *_):
        obs = adata.obs
        m = np.ones(adata.n_obs, bool)
        if self.f_donor.value != self.ALL:
            m &= (obs["donor"] == self.f_donor.value).values
        if self.f_ct.value != self.ALL:
            m &= (obs["cell_type"] == self.f_ct.value).values
        if self.f_cyt.value != self.ALL:
            m &= (obs["cytokine"] == self.f_cyt.value).values
        matches = [s for s, p in zip(self.ids, self.pos) if m[p]]
        shown = matches[:self.max_options]
        opts = sorted(set(shown) | set(self.chosen), key=self.rank.get)
        self._busy = True
        try:
            self.select.value = ()
            self.select.options = opts
            self.select.value = tuple(s for s in self.chosen if s in opts)
        finally:
            self._busy = False
        extra = f" (showing the first {self.max_options}: narrow it down with the filters)" if len(matches) > self.max_options else ""
        self._info = f"{len(matches):,} samples match the filters{extra}. Ctrl/⌘ + click selects several."
        self._update_status()

    def _update_status(self):
        self.status.value = f'<div class="cx-hint">{self._info} <b>{len(self.chosen)}</b> selected.</div>'

    def _on_select(self, change):
        if self._busy:
            return
        self.chosen = list(change["new"])
        self.value_widget.value = "\n".join(self.chosen)
        self._update_status()


def gene_picker(default, width="330px", placeholder="STAT1, CXCL9, ...", with_sets=True, label="Genes (gene symbols, separated by commas)"):
    """Genes as a list of symbols (editable text) + a ready-made-lists dropdown + an auto-complete box to add one gene."""
    text = text_input(default, width, placeholder)
    sets = geneset_dropdown(text) if with_sets else None
    add = widgets.Combobox(value="", options=list(GENE_NAMES), ensure_option=False, continuous_update=False,
                           placeholder="type a gene, press Enter", layout=_lay("190px"), style=STYLE)

    def _add(change):
        new = change["new"]
        if new:
            g, _unknown = parse_genes(new, 1)
            if g and g[0] not in parse_genes(text.value, 500)[0]:
                text.value = (text.value.strip().rstrip(",") + ", " + g[0]).lstrip(", ")
            add.value = ""
    add.observe(_add, "value")
    second = [field("…or a ready-made list", sets)] if with_sets else []
    return text, widgets.VBox([field(label, text), widgets.HBox(second + [field("…or add one gene" if with_sets else "Add one gene (auto-complete)", add)])])


MATRICES = [("X: raw counts", "X"), ("layer 'cpm': counts per million", "cpm"), ("layer 'log2cpm': log2(CPM + 1)", "log2cpm")]


def anndata_repr(a):
    """The classic text of an AnnData object (we write it ourselves so that it looks the same with every anndata version)."""
    lines = [f"AnnData object with n_obs × n_vars = {a.n_obs} × {a.n_vars}",
             "    obs: " + ", ".join(f"'{c}'" for c in a.obs.columns),
             "    var: " + ", ".join(f"'{c}'" for c in a.var.columns)]
    layers = [k for k in a.layers.keys() if k is not None]
    if layers:
        lines.append("    layers: " + ", ".join(f"'{k}'" for k in layers))
    return "\n".join(lines)


def _anndata_cartoon(sub):
    """Left: the AnnData scheme with the shapes of the selection. Right: where the selection sits in the full matrix."""
    n, m = sub.n_obs, sub.n_vars
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.3), gridspec_kw={"width_ratios": [1.55, 1]})
    a1.set_xlim(0, 10); a1.set_ylim(0, 5.6); a1.axis("off")
    box = lambda x, y, w, h, col, txt, tc="white": (a1.add_patch(plt.Rectangle((x, y), w, h, fc=col, ec="white", lw=2)),
                                                     a1.text(x + w / 2, y + h / 2, txt, ha="center", va="center", color=tc, fontsize=10, fontweight="bold"))
    box(2.5, 0.4, 5.2, 3.5, BLUE, f"X\n{n} samples × {m} genes\n(the numbers)")
    box(0.1, 0.4, 2.2, 3.5, ORANGE, f"obs\n{n} samples ×\n{sub.obs.shape[1]} columns\n(about the\nSAMPLES)")
    box(2.5, 4.1, 5.2, 1.3, GREEN, f"var\n{m} genes × {sub.var.shape[1]} columns\n(about the GENES)")
    a1.text(8.0, 2.3, f"layers\n{', '.join(k for k in sub.layers.keys() if k is not None)}\n(same shape as X)", va="center", fontsize=8.5, color=PBS_COLOR)
    a1.set_title("The AnnData object you selected", loc="left")
    # right: where the subset is in the complete matrix
    a2.set_xlim(0, adata.n_vars); a2.set_ylim(adata.n_obs, 0)
    a2.add_patch(plt.Rectangle((0, 0), adata.n_vars, adata.n_obs, fc="#eef3fb", ec=BLUE, lw=1.2))
    ri, ci = adata.obs_names.get_indexer(sub.obs_names), adata.var_names.get_indexer(sub.var_names)
    for r in ri:
        a2.axhline(r, color=ORANGE, lw=1.4, alpha=0.9)
    for c in ci:
        a2.axvline(c, color=GREEN, lw=1.4, alpha=0.9)
    a2.scatter(np.tile(ci, len(ri)), np.repeat(ri, len(ci)), s=14, color=BLUE, zorder=3)
    a2.set_xticks([]); a2.set_yticks([])
    a2.set_xlabel(f"{adata.n_vars:,} genes", fontsize=9); a2.set_ylabel(f"{adata.n_obs:,} samples", fontsize=9)
    a2.spines[:].set_visible(False)
    a2.set_title("…inside the complete matrix", loc="left")
    fig.tight_layout()
    show_fig(fig)


def ex3_render(samples, genes_text, matrix):
    genes, unknown = parse_genes(genes_text, 200)
    unknown_genes_message(unknown)
    samples = [s for s in str(samples).split("\n") if s in adata.obs_names]
    if not samples or not genes:
        msg("Select <b>at least one sample</b> and <b>one gene</b> to create a subset.", "info")
        return
    sub = adata[samples, genes]                                                    # ← the subset: samples first, genes second
    display(HTML(f'<div class="cx-label">The code that creates your subset</div>'
                 f'<div class="cx-code">sub = adata[samples, genes]      <span class="cx-com"># {len(samples)} sample names, {len(genes)} gene names</span>\n'
                 f'print(sub)</div>'
                 f'<div class="cx-code" style="margin-top:4px;background:rgba(42,120,214,.10)">{html.escape(anndata_repr(sub))}</div>'))
    _anndata_cartoon(sub)
    M = sub.X if matrix == "X" else sub.layers[matrix]
    dfx = pd.DataFrame(np.asarray(M), index=sub.obs_names, columns=sub.var_names)
    fmt = "{:,.0f}" if matrix == "X" else ("{:,.1f}" if matrix == "cpm" else "{:.2f}")
    shade = np.log10(dfx.values.astype(float) + 1) if matrix != "log2cpm" else dfx.values.astype(float)
    top = max(float(np.nanmax(shade)), 1e-9) * 1.25
    head = "".join(f'<th>{html.escape(str(c))}</th>' for c in dfx.columns)                      # plain HTML: no jinja2 needed
    body = "".join(
        f'<tr><td class="l"><b>{html.escape(str(idx))}</b></td>' + "".join(
            f'<td style="background:rgba(42,120,214,{0.05 + 0.6 * max(0.0, float(sh)) / top:.2f})">{fmt.format(v)}</td>' for v, sh in zip(row, shrow)) + "</tr>"
        for idx, row, shrow in zip(dfx.index, dfx.values, shade))
    table_x = f'<table class="cx-table"><thead><tr><th class="l"></th>{head}</tr></thead><tbody>{body}</tbody></table>'
    obs_df = sub.obs.copy()
    for c in ("donor", "cytokine", "treatment", "cell_type"):
        obs_df[c] = obs_df[c].astype(str)
    obs_df["psbulk_n_cells"] = obs_df["psbulk_n_cells"].map("{:,}".format); obs_df["library_size"] = obs_df["library_size"].map("{:,}".format)
    var_df = sub.var.copy()
    var_df["total_counts"] = var_df["total_counts"].map("{:,}".format)
    tabs = widgets.Output(), widgets.Output(), widgets.Output()
    with tabs[0]:
        display(HTML(f'<div class="cx-hint">Rows = the {len(samples)} samples (same rows as <code>obs</code>), columns = the {len(genes)} genes (same rows as <code>var</code>).</div>'))
        display(HTML(table_x))
    with tabs[1]:
        display(HTML('<div class="cx-hint">One row per sample: who it comes from, what was added, which cell type, and how many cells / counts it holds.</div>'))
        display(HTML(html_table(obs_df, left_cols=obs_df.columns, index_label="sample ID")))
    with tabs[2]:
        display(HTML('<div class="cx-hint">One row per gene: its total counts in the complete dataset, and in how many of the 4,224 samples it has at least one count.</div>'))
        display(HTML(html_table(var_df, left_cols=["total_counts"], index_label=disp("gene"))))
    tab = widgets.Tab(children=list(tabs))
    for i, t in enumerate(["X  (the numbers)", "obs  (samples)", "var  (genes)"]):
        tab.set_title(i, t)
    display(tab)


def exercise_anndata():
    default = [f"Donor1|{ct}|PBS" for ct in CELL_TYPES] + ["Donor1|CD14_Mono|IFN-gamma"]
    sp = SamplePicker(default)
    gtext, gwidget = gene_picker("LYZ, CD3E, NKG7, STAT1, ACTB")
    matrix = dropdown(MATRICES, "X", "300px")
    c = dict(samples=sp.value_widget, genes_text=gtext, matrix=matrix)
    rows = [widgets.HBox([widgets.VBox([widgets.HTML('<div class="cx-label">① Choose samples (rows)</div>'), sp.widget]),
                          widgets.VBox([widgets.HTML('<div class="cx-label">② Choose genes (columns)</div>'), gwidget,
                                        field("③ Which matrix to look at", matrix, "X holds the raw counts; the layers hold normalised versions (Exercises 5-6)")],
                                       layout=_lay(margin="0 0 0 20px"))],
                         layout=_lay(flex_flow="row wrap"))]
    ex = Exercise(c, ex3_render, rows, title="Build your own subset of the dataset")
    ex.picker = sp
    ex.show()
    return ex
