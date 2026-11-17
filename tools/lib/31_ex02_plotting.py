# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  Exercise 2: basic plots (simulated data, independent of the dataset)
# ══════════════════════════════════════════════════════════════════════════════════════════════════
PLOT_TYPES = [("📊 Histogram", "hist"), ("📦 Box plot", "box"), ("📶 Bar chart", "bar"), ("🔵 Scatter plot", "scatter"), ("📈 Line plot (time course)", "line")]
COL_CTRL, COL_TRT = "#8a8984", "#eb6834"


def sim_expression(n, mean, effect, sd, r, seed):
    """Fake expression values (log2 CPM) of ONE gene in `n` donors: without (control) and with (treated) a stimulus.
    `effect` = log2 fold change; `r` = how much a donor's two values resemble each other (same person, same baseline)."""
    rng = np.random.default_rng(int(seed))
    z1 = rng.standard_normal(n)
    z2 = r * z1 + np.sqrt(max(1 - r ** 2, 0)) * rng.standard_normal(n)
    return mean + sd * z1, mean + effect + sd * z2


def ex2_render(plot_type, n, mean, effect, sd, bins, r, tau, seed):
    control, treated = sim_expression(n, mean, effect, sd, r, seed)
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    if plot_type == "hist":
        lo, hi = min(control.min(), treated.min()), max(control.max(), treated.max())
        edges = np.linspace(lo, hi, bins + 1)
        ax.hist(control, bins=edges, color=COL_CTRL, alpha=0.65, label=disp("control"))
        ax.hist(treated, bins=edges, color=COL_TRT, alpha=0.65, label=disp("treated"))
        ax.axvline(control.mean(), color=COL_CTRL, lw=2, ls="--"); ax.axvline(treated.mean(), color=COL_TRT, lw=2, ls="--")
        ax.set_xlabel("expression (log2 CPM)"); ax.set_ylabel("number of donors"); ax.legend()
        ax.set_title(f"Histogram: how many donors have each value? (dashed = mean)")
        code = f"plt.hist(treated, bins={bins})"
    elif plot_type == "box":
        bp = ax.boxplot([control, treated], widths=0.5, patch_artist=True, showfliers=False, medianprops=dict(color="k", lw=2))
        for patch, col in zip(bp["boxes"], (COL_CTRL, COL_TRT)):
            patch.set(facecolor=col, alpha=0.45)
        rng = np.random.default_rng(1)
        for i, v in enumerate((control, treated), 1):
            ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=14, color=(COL_CTRL, COL_TRT)[i - 1], alpha=0.8, zorder=3)
        ax.set_xticks([1, 2]); ax.set_xticklabels([disp("control"), disp("treated")]); ax.set_ylabel("expression (log2 CPM)")
        ax.set_title("Box plot: median (line), middle 50 % (box), range (whiskers)")
        code = "plt.boxplot([control, treated])"
    elif plot_type == "bar":
        means, sds = [control.mean(), treated.mean()], [control.std(ddof=1), treated.std(ddof=1)]
        ax.bar([0, 1], means, yerr=sds, capsize=8, color=(COL_CTRL, COL_TRT), alpha=0.8, width=0.55)
        ax.set_xticks([0, 1]); ax.set_xticklabels([disp("control"), disp("treated")]); ax.set_ylabel("expression (log2 CPM): mean ± SD")
        ax.set_title("Bar chart: one number (the mean) per group")
        code = "plt.bar(['control', 'treated'], [control.mean(), treated.mean()])"
    elif plot_type == "scatter":
        ax.scatter(control, treated, s=26, color=BLUE, alpha=0.75, edgecolor="white", linewidth=0.5)
        lo, hi = min(control.min(), treated.min()) - 0.5, max(control.max(), treated.max()) + 0.5
        ax.plot([lo, hi], [lo, hi], color="k", lw=1, ls="--"); ax.text(hi, hi, " y = x", va="top", ha="right", fontsize=8)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
        rr = np.corrcoef(control, treated)[0, 1] if n > 2 else float("nan")
        ax.set_xlabel("control (log2 CPM)"); ax.set_ylabel("treated (log2 CPM)")
        ax.set_title(f"Scatter plot: one point per donor (correlation r = {rr:.2f})")
        code = "plt.scatter(control, treated)"
    else:                                                          # line plot: a time course
        t = np.linspace(0, 24, 13)
        rng = np.random.default_rng(int(seed) + 1000)
        k = min(int(n), 12)                                         # at most 12 donors as thin lines
        for i in range(k):
            ax.plot(t, mean + (control[i] - mean) + effect * (1 - np.exp(-t / tau)) + rng.normal(0, 0.12 * sd, len(t)), color=COL_TRT, alpha=0.25, lw=1)
        ax.plot(t, mean + effect * (1 - np.exp(-t / tau)), color=COL_TRT, lw=3, label="mean response (treated)")
        ax.axhline(mean, color=COL_CTRL, lw=2, ls="--", label="control (no change)")
        ax.set_xlabel("hours after adding the cytokine"); ax.set_ylabel("expression (log2 CPM)"); ax.legend(loc="best")
        ax.set_title(f"Line plot: expression over time ({k} donors, thin lines)")
        code = "plt.plot(hours, expression)"
    fig.tight_layout()
    show_fig(fig)
    display(HTML(f'<div class="cx-hint">Behind the plot, Python runs just one line:</div><div class="cx-code">{html.escape(code)}</div>'))


def exercise_plotting():
    c = dict(plot_type=dropdown(PLOT_TYPES, "hist", "230px"), n=int_slider(5, 200, 60), mean=float_slider(0, 10, 5, 0.5, fmt=".1f"),
             effect=float_slider(-3, 3, 1.0, 0.25, fmt="+.2f"), sd=float_slider(0.1, 3, 1.0, 0.1, fmt=".1f"), bins=int_slider(3, 60, 20),
             r=float_slider(0, 0.99, 0.7, 0.05), tau=float_slider(1, 12, 4, 1, fmt=".0f"), seed=widgets.IntText(value=1, layout=_lay("70px")))
    dice = widgets.Button(description="🎲 new random draw", layout=_lay("170px"))
    dice.on_click(lambda _b: setattr(c["seed"], "value", c["seed"].value + 1))
    f_bins = field("bins  (histogram)", c["bins"], "how many bars")
    f_r = field("similarity of the two values of a donor  (scatter)", c["r"], "0 = unrelated, 1 = identical pattern")
    f_tau = field("response speed: hours to reach 63 %  (line)", c["tau"])
    rows = [[field("type of plot", c["plot_type"]), field("number of donors", c["n"]), field("mean expression (control)", c["mean"])],
            [field("effect of the cytokine (log2 fold change)", c["effect"], "+1 = twice as much, −1 = half"), field("variability between donors (SD)", c["sd"]),
             f_bins, f_r, f_tau],
            [field("random draw number", c["seed"]), widgets.VBox([widgets.HTML('<div class="cx-label">&nbsp;</div>'), dice])]]
    show_for = {"hist": [f_bins], "scatter": [f_r], "line": [f_tau]}

    def _toggle(*_):
        for f in (f_bins, f_r, f_tau):
            f.layout.display = "none"
        for f in show_for.get(c["plot_type"].value, []):
            f.layout.display = None
    c["plot_type"].observe(_toggle, "value")
    _toggle()
    ex = Exercise(c, ex2_render, rows, title="Simulate a gene and choose how to look at it",
                  note="The data are <b>simulated</b> (no real dataset yet): a gene measured in several donors, without and with a cytokine.")
    ex.show()
    return ex
