# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  20_widgets.py: the small toolbox every exercise is built from (controls, validation, layout, figures)
# ══════════════════════════════════════════════════════════════════════════════════════════════════
STRICT = False            # tests set this to True so that errors are raised instead of shown as a friendly message
FIG_SINK = None           # tests set this to a function(fig) that receives every figure shown (to save a PNG)
ERRORS = []               # every error caught while redrawing (the students see a friendly message; the tests read this list)

CSS = """
<style>
.cx-panel {border:1px solid rgba(127,127,127,.40); border-radius:8px; padding:8px 12px 6px 12px; background:rgba(127,127,127,.08); margin:2px 0 8px 0;}
.cx-title {font-weight:700; font-size:14px; margin:0 0 4px 0;}
.cx-label {font-weight:600; font-size:12px; margin:0 0 1px 2px;}
.cx-hint  {font-size:11px; opacity:.7; margin:0 0 2px 2px;}
.cx-warn  {border-left:4px solid #eda100; background:rgba(237,161,0,.14); padding:6px 10px; margin:4px 0; border-radius:4px; font-size:13px;}
.cx-err   {border-left:4px solid #e34948; background:rgba(227,73,72,.12); padding:6px 10px; margin:4px 0; border-radius:4px; font-size:13px;}
.cx-info  {border-left:4px solid #2a78d6; background:rgba(42,120,214,.10); padding:6px 10px; margin:4px 0; border-radius:4px; font-size:13px;}
.cx-code  {font-family:monospace; font-size:13px; line-height:1.45; white-space:pre; border:1px solid rgba(127,127,127,.35); border-radius:6px; padding:8px 12px; background:rgba(127,127,127,.07); overflow-x:auto;}
.cx-kw {color:#7a4fd6; font-weight:600} .cx-str {color:#1a8f5a} .cx-num {color:#d6571f} .cx-com {color:#8a8984; font-style:italic}
table.cx-table {border-collapse:collapse; font-size:12.5px; margin:2px 0;}
table.cx-table th {text-align:right; padding:3px 9px; border-bottom:2px solid rgba(127,127,127,.5); font-weight:600;}
table.cx-table td {text-align:right; padding:2px 9px; border-bottom:1px solid rgba(127,127,127,.2);}
table.cx-table td.l, table.cx-table th.l {text-align:left;}
table.cx-table tr:hover td {background:rgba(127,127,127,.12);}
.cx-up {color:#e34948; font-weight:600} .cx-down {color:#2a78d6; font-weight:600}
</style>
"""
STYLE = {"description_width": "initial"}


# ── messages shown in the output area (never a traceback) ───────────────────────────────────────────
def msg(text, kind="warn"):
    """Show a coloured message box. kind: 'warn' (orange), 'err' (red) or 'info' (blue)."""
    display(HTML(f'<div class="cx-{kind}">{text}</div>'))


# ── control builders ───────────────────────────────────────────────────────────────────────────────
def _lay(width=None, **kw):
    return widgets.Layout(width=width, **kw) if width else widgets.Layout(**kw)


def field(label, widget, hint=None):
    """A control with a small bold label above it (and an optional grey hint below)."""
    parts = [widgets.HTML(f'<div class="cx-label">{label}</div>'), widget]
    if hint:
        parts.append(widgets.HTML(f'<div class="cx-hint">{hint}</div>'))
    return widgets.VBox(parts, layout=_lay(margin="0 14px 4px 0"))


def dropdown(options, value=None, width="230px"):
    return widgets.Dropdown(options=options, value=value, layout=_lay(width), style=STYLE)


def toggle(options, value=None, width=None):
    return widgets.ToggleButtons(options=options, value=value, style={"button_width": width or "auto"})


def int_slider(lo, hi, value, step=1, width="260px"):
    return widgets.IntSlider(min=lo, max=hi, value=value, step=step, continuous_update=False, layout=_lay(width), style=STYLE)


def float_slider(lo, hi, value, step=0.1, width="260px", fmt=".2f"):
    return widgets.FloatSlider(min=lo, max=hi, value=value, step=step, continuous_update=False, readout_format=fmt, layout=_lay(width), style=STYLE)


def checkbox(label, value=False, width="260px"):
    return widgets.Checkbox(value=value, description=label, indent=False, layout=_lay(width))


def text_input(value="", width="320px", placeholder=""):
    return widgets.Text(value=value, placeholder=placeholder, continuous_update=False, layout=_lay(width), style=STYLE)


def multi_select(options, value, rows=6, width="200px"):
    return widgets.SelectMultiple(options=options, value=tuple(value), rows=rows, layout=_lay(width), style=STYLE)


def selection_summary(select, label_fn=lambda v: v, max_items=8):
    """A small line under a multi-select that spells out what is selected (the list only shows its first rows)."""
    box = widgets.HTML()

    def _update(*_):
        vals = list(select.value)
        shown = ", ".join(html.escape(str(label_fn(v))) for v in vals[:max_items]) + (f" … (+{len(vals) - max_items})" if len(vals) > max_items else "")
        box.value = f'<div class="cx-hint"><b>{len(vals)}</b> selected: {shown if vals else "none"}</div>'
    select.observe(_update, "value")
    _update()
    return box


def gene_combo(value="ACTB", width="200px"):
    """One gene symbol with auto-completion (type a few letters, pick from the list, press Enter)."""
    return widgets.Combobox(value=value, options=list(GENE_NAMES), ensure_option=False, continuous_update=False,
                            placeholder="gene symbol", layout=_lay(width), style=STYLE)


def geneset_dropdown(text_widget, default_label="— choose a ready-made list —", width="260px"):
    """A dropdown of gene sets that, when used, writes the genes of the set into `text_widget`."""
    dd = dropdown([(default_label, "")] + [(k, k) for k in GENESETS], "", width)

    def _on(change):
        if change["new"]:
            text_widget.value = ", ".join(GENESETS[change["new"]])
            dd.value = ""                                         # back to the placeholder, so the same set can be chosen again
    dd.observe(_on, "value")
    return dd


def panel(title, *rows, note=None):
    """Titled grey box with rows of controls (each row = a widget or a list of widgets)."""
    kids = [widgets.HTML(f'<div class="cx-title">🎛 {title}</div>')]
    for r in rows:
        kids.append(widgets.HBox(list(r), layout=_lay(flex_flow="row wrap", align_items="flex-start")) if isinstance(r, (list, tuple)) else r)
    if note:
        kids.append(widgets.HTML(f'<div class="cx-hint">{note}</div>'))
    box = widgets.VBox(kids)
    box.add_class("cx-panel")
    return box


class Exercise:
    """Controls + a function that redraws the outputs every time a control changes (`widgets.interactive_output`).
    `render(**values)` receives the current value of every control (by name) and draws with matplotlib / display()."""

    def __init__(self, controls, render, ui_rows, title="Controls", note=None, extra=()):
        self.render, self.paused = render, False
        self.controls = dict(controls)
        self._tick = widgets.IntText(value=0)                         # (never displayed) bumping it redraws once, after several changes
        self.panel = panel(title, *ui_rows, note=note)
        self.output = widgets.interactive_output(self._call, {**self.controls, "_tick": self._tick})
        self.ui = widgets.VBox([widgets.HTML(CSS), self.panel, *extra, self.output])

    def values(self):
        return {k: w.value for k, w in self.controls.items()}

    def set_values(self, **values):
        """Change several controls and redraw ONCE."""
        self.paused = True
        try:
            for k, v in values.items():
                self.controls[k].value = v
        finally:
            self.paused = False
        self._tick.value += 1

    def _call(self, _tick=0, **kw):
        if self.paused:
            return
        try:
            self.render(**kw)
        except Exception as err:                                  # a student must never see a traceback
            ERRORS.append(f"{type(err).__name__}: {err}")
            if STRICT:
                raise
            msg(f"⚠️ Something went wrong with this combination ({html.escape(type(err).__name__)}: {html.escape(str(err))[:200]}). Try changing a control, or go back to the default values.", "err")

    def show(self):
        display(self.ui)
        return None


# ── figures ────────────────────────────────────────────────────────────────────────────────────────
def show_fig(fig):
    """Show a figure in the current output and free its memory."""
    if FIG_SINK is not None:
        FIG_SINK(fig)
    plt.show()
    plt.close(fig)


def fmt_n(x):
    return f"{x:,.0f}"


def si(x):
    """12,345,678 → '12.3 M'"""
    x = float(x)
    for lim, s in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if abs(x) >= lim:
            return f"{x / lim:.1f} {s}"
    return f"{x:.0f}"


def html_table(df, index=True, float_fmt="{:,.2f}", classes="cx-table", left_cols=(), index_label=None, max_rows=200):
    """Small, clean HTML table (works in Jupyter, Colab, VS Code)."""
    df = df.head(max_rows)
    cols = list(df.columns)
    head = "".join(f'<th class="{"l" if c in left_cols else ""}">{html.escape(str(c))}</th>' for c in cols)
    if index:
        head = f'<th class="l">{html.escape(index_label or "")}</th>' + head
    rows = []
    for idx, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (float, np.floating)):
                s = "—" if not np.isfinite(v) else (float_fmt.format(v))
            elif isinstance(v, (int, np.integer)):
                s = f"{int(v):,}"
            else:
                s = str(v) if str(v).startswith("<") else html.escape(str(v))          # allow our own <span> markup
            cells.append(f'<td class="{"l" if c in left_cols else ""}">{s}</td>')
        rows.append("<tr>" + (f'<td class="l"><b>{html.escape(str(idx))}</b></td>' if index else "") + "".join(cells) + "</tr>")
    return f'<table class="{classes}"><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table>'


def show_table(df, **kw):
    display(HTML(html_table(df, **kw)))


# ── genes: parsing and validation (names are case-insensitive; unknown ones get a friendly hint) ─────────────
GENE_NAMES, GENE_LOOKUP, GENESETS = [], {}, {}


def parse_genes(text, max_genes=60):
    """'stat1, Ifit3  nonsense' → (['STAT1', 'IFIT3'], ['nonsense'])"""
    valid, unknown = [], []
    for tok in [t for t in re.split(r"[,\s;]+", str(text).strip()) if t]:
        name = GENE_LOOKUP.get(tok.upper())
        if name is None:
            unknown.append(tok)
        elif name not in valid:
            valid.append(name)
    return valid[:max_genes], unknown


def unknown_genes_message(unknown):
    if not unknown:
        return
    parts = []
    for t in unknown[:6]:
        close = difflib.get_close_matches(t.upper(), list(GENE_LOOKUP), n=3, cutoff=0.7)
        parts.append(f"<b>{html.escape(t)}</b>" + (f" (did you mean {', '.join(GENE_LOOKUP[c] for c in close)}?)" if close else ""))
    msg("Gene not found in this dataset: " + "; ".join(parts) + ". Check the spelling (gene symbols, e.g. <i>STAT1</i>).", "warn")


def one_gene(text, default=None):
    """Validate a single gene symbol; shows a hint and returns `default` (possibly None) if it is not valid."""
    valid, unknown = parse_genes(text, 1)
    unknown_genes_message(unknown if not valid else [])
    return valid[0] if valid else default


def gene_index(gene):
    return adata.var_names.get_loc(gene)


def sample_mask(cell_types=None, cytokines=None, donors=None):
    """Boolean mask over the samples of `adata` (None = no restriction)."""
    m = np.ones(adata.n_obs, bool)
    if cell_types is not None:
        m &= adata.obs["cell_type"].isin(list(cell_types)).values
    if cytokines is not None:
        m &= adata.obs["cytokine"].isin(list(cytokines)).values
    if donors is not None:
        m &= adata.obs["donor"].isin(list(donors)).values
    return m


def cyt_label(c):
    """Cytokine names as they are written in the data (IFN-gamma, IL-1-beta) → nicer for titles: IFN-γ, IL-1β."""
    c = re.sub(r"(IL-\d+)-(alpha|beta)", lambda m: m.group(1) + {"alpha": "α", "beta": "β"}[m.group(2)], str(c))        # IL-1-beta → IL-1β
    return c.replace("-gamma", "-γ").replace("-beta", "-β").replace("-alpha", "-α").replace("-omega", "-ω").replace("-lambda", "-λ").replace("-epsilon", "-ε")


def pick_default(options, wanted, fallback_n=3):
    """The wanted values that exist in `options` (so a default never breaks if the data change)."""
    ok = [w for w in wanted if w in options]
    return ok or list(options)[:fallback_n]
