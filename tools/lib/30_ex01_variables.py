# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  Exercise 1: variables and print()
# ══════════════════════════════════════════════════════════════════════════════════════════════════
def highlight_python(code):
    """Tiny syntax highlighter (strings, numbers, comments, keywords) → HTML for the code panel."""
    out = []
    pattern = re.compile(r"(#.*$)|(f?'[^']*'|f?\"[^\"]*\")|(\b\d+(?:\.\d+)?\b)|(\b(?:print|len|True|False)\b)")
    for line in html.escape(code, quote=False).split("\n"):
        def sub(m):
            for cls, g in zip(("com", "str", "num", "kw"), m.groups()):
                if g:
                    return f'<span class="cx-{cls}">{g}</span>'
        out.append(pattern.sub(sub, line))
    return "\n".join(out)


def ex1_program(message, n_donors, cytokine, genes_text, is_control, n_cells, umis_per_cell):
    """The Python program that corresponds to the current values of the controls (as text)."""
    genes = [g.strip() for g in re.split(r"[,;]", genes_text) if g.strip()]
    return "\n".join([
        "# 1) Create variables: a NAME on the left of  =  and a VALUE on the right",
        f"message = {message!r}",
        f"n_donors = {int(n_donors)}",
        f"cytokine = {cytokine!r}",
        f"genes = {genes!r}",
        f"is_control = {bool(is_control)}",
        f"n_cells = {int(n_cells)}",
        f"umis_per_cell = {float(umis_per_cell)!r}",
        "",
        "# 2) Use them: print() shows a value on the screen",
        "print(message)",
        "print('Donors:', n_donors)",
        "print(cytokine, 'was added to', n_donors, 'donors')",
        "print('Number of genes:', len(genes))",
        "",
        "# 3) Variables can be combined to create new variables",
        "library_size = n_cells * umis_per_cell",
        "print('Library size:', library_size)",
        "print('Is this sample a control?', is_control)",
    ])


def ex1_run(program):
    """Run the (safe: built by ex1_program from repr() literals) program; returns what it printed and its variables."""
    buf, ns = io.StringIO(), {}
    real_stdout = sys.stdout
    sys.stdout = buf
    try:
        exec(program, {"__builtins__": {"print": print, "len": len}}, ns)
    finally:
        sys.stdout = real_stdout
    return buf.getvalue(), ns


def ex1_render(message, n_donors, cytokine, genes_text, is_control, n_cells, umis_per_cell):
    program = ex1_program(message, n_donors, cytokine, genes_text, is_control, n_cells, umis_per_cell)
    printed, ns = ex1_run(program)
    badge = {"str": "#1a8f5a", "int": "#d6571f", "float": "#d6571f", "bool": "#7a4fd6", "list": "#2a78d6"}
    rows = "".join(
        f'<tr><td class="l"><b>{k}</b></td><td class="l">{html.escape(repr(v))}</td>'
        f'<td class="l"><span style="color:{badge.get(type(v).__name__, "#888")};font-weight:600">{type(v).__name__}</span></td></tr>'
        for k, v in ns.items())
    table = f'<table class="cx-table"><thead><tr><th class="l">name</th><th class="l">value stored</th><th class="l">type</th></tr></thead><tbody>{rows}</tbody></table>'
    display(HTML(
        '<div style="display:flex; flex-wrap:wrap; gap:16px; align-items:flex-start">'
        f'<div style="flex:1 1 440px"><div class="cx-label">🐍 The Python code written by your controls</div><div class="cx-code">{highlight_python(program)}</div></div>'
        '<div style="flex:1 1 320px">'
        f'<div class="cx-label">🖨️ What <code>print()</code> shows</div><div class="cx-code" style="background:rgba(42,120,214,.10)">{html.escape(printed)}</div>'
        f'<div class="cx-label" style="margin-top:10px">🗃️ What Python remembers (its variables)</div>{table}</div></div>'))


def exercise_variables():
    c = dict(
        message=text_input("Hello, immunologists!", "300px"),
        n_donors=int_slider(1, 12, 12, width="240px"),
        cytokine=dropdown(["IFN-gamma", "IL-4", "IL-6", "TNF-alpha", "IL-2", "IL-10"], "IFN-gamma", "170px"),
        genes_text=text_input("STAT1, CXCL9, ACTB", "250px"),
        is_control=checkbox("this sample is a control", False, "210px"),
        n_cells=int_slider(100, 5000, 2500, 100, "240px"),
        umis_per_cell=float_slider(1000, 8000, 3000, 250, "240px", ".0f"),
    )
    rows = [[field("message  (text)", c["message"], "any text"), field("cytokine  (text)", c["cytokine"]),
             field("genes  (a list of texts)", c["genes_text"], "separate with commas")],
            [field("n_donors  (whole number)", c["n_donors"]), field("n_cells  (whole number)", c["n_cells"]),
             field("umis_per_cell  (decimal number)", c["umis_per_cell"])],
            [field("is_control  (True / False)", c["is_control"])]]
    ex = Exercise(c, ex1_render, rows, title="Create your variables")
    ex.show()
    return ex
