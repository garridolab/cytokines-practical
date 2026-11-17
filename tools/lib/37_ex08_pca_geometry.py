# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  Exercise 8: the geometry of PCA with two genes (rotation and projection)
# ══════════════════════════════════════════════════════════════════════════════════════════════════
PAIRS8 = {"corr": ("ISG15", "MX1"), "uncorr": ("IDO1", "IL1B")}
PAIR_CHOICES = [("Correlated genes: ISG15 & MX1 (two interferon-response genes)", "corr"), ("Uncorrelated genes: IDO1 & IL1B (two different pathways)", "uncorr")]


def pair_data(pair):
    """Centred log2(CPM + 1) of the two genes in the monocyte samples + everything PCA would compute from them."""
    a, b = PAIRS8[pair]
    m = adata.obs["cell_type"].values == "CD14_Mono"
    P = np.c_[adata.layers["log2cpm"][m, gene_index(a)], adata.layers["log2cpm"][m, gene_index(b)]].astype(np.float64)
    Pc = P - P.mean(axis=0)
    C = np.cov(Pc.T)
    w, V = np.linalg.eigh(C)                                                  # PCA of two variables: eigen-decomposition of the covariance
    w, V = w[::-1], V[:, ::-1]
    if V[np.argmax(np.abs(V[:, 0])), 0] < 0:
        V[:, 0] = -V[:, 0]
    theta1 = np.degrees(np.arctan2(V[1, 0], V[0, 0])) % 180
    is_ifn = adata.obs.loc[adata.obs["cell_type"] == "CD14_Mono", "cytokine"].astype(str).str.startswith("IFN-").values
    return dict(genes=(a, b), Pc=Pc, C=C, eigvals=w, V=V, theta1=theta1, r=np.corrcoef(Pc.T)[0, 1], is_ifn=is_ifn)


def ex8_render(pair, angle, show_proj, show_pc):
    d = pair_data(pair)
    a, b = d["genes"]
    Pc, C = d["Pc"], d["C"]
    th = np.radians(angle)
    u = np.array([np.cos(th), np.sin(th)])                                   # direction of our axis
    v = np.array([-np.sin(th), np.cos(th)])                                  # the perpendicular one
    s_u, s_v = Pc @ u, Pc @ v                                                # coordinates in the rotated system
    total = np.trace(C)
    captured = u @ C @ u
    best = d["eigvals"][0]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), gridspec_kw={"width_ratios": [1.15, 1, 1.15]})
    # ① the data with the axis and the projections
    ax = axes[0]
    ax.scatter(Pc[:, 0], Pc[:, 1], s=7, alpha=0.4, color=NS_COLOR, linewidths=0, label="other cytokines and PBS")
    ax.scatter(Pc[d["is_ifn"], 0], Pc[d["is_ifn"], 1], s=9, alpha=0.7, color=UP_COLOR, linewidths=0, label="interferons")
    L = 1.15 * np.abs(Pc).max()
    if show_proj:
        idx = np.random.default_rng(3).choice(len(Pc), 45, replace=False)
        for i in idx:
            p = (Pc[i] @ u) * u
            ax.plot([Pc[i, 0], p[0]], [Pc[i, 1], p[1]], color="k", lw=0.7, ls=":", alpha=0.8)
        ax.scatter(*(np.outer(Pc[idx] @ u, u)).T, s=12, color="k", zorder=5, label="projection on the axis")
    ax.plot([-L * u[0], L * u[0]], [-L * u[1], L * u[1]], color=PURPLE, lw=2.6, zorder=4, label=f"your axis ({angle}°)")
    if show_pc:
        for k, (col, nm) in enumerate(((DARK, "PC1"), ("#888", "PC2"))):
            e = d["V"][:, k] * 1.9 * np.sqrt(d["eigvals"][k])
            ax.annotate("", xy=e, xytext=-e, arrowprops=dict(arrowstyle="<|-|>", color=col, lw=1.6, ls="--"))
            ax.text(*(e * 1.08), nm, fontsize=10, fontweight="bold", color=col, ha="center", va="center")
    ax.axhline(0, color="#ccc", lw=0.6); ax.axvline(0, color="#ccc", lw=0.6)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal")
    ax.set_xlabel(f"{a}  (log₂ CPM + 1, centred)"); ax.set_ylabel(f"{b}  (log₂ CPM + 1, centred)")
    ax.set_title(f"① {a} vs {b} in monocytes (r = {d['r']:.2f})")
    ax.legend(loc="lower right", fontsize=7.5, markerscale=1.4)
    # ② variance of the projections for every possible angle
    ax = axes[1]
    ang = np.arange(0, 181)
    var_ang = np.array([np.array([np.cos(np.radians(t)), np.sin(np.radians(t))]) @ C @ np.array([np.cos(np.radians(t)), np.sin(np.radians(t))]) for t in ang])
    ax.plot(ang, var_ang / total * 100, color=BLUE, lw=2.5)
    ax.scatter([angle], [captured / total * 100], s=90, color=PURPLE, zorder=5)
    ax.annotate(f"{captured / total:.0%}", xy=(angle, captured / total * 100), xytext=(0, 11), textcoords="offset points", ha="center", fontsize=10.5, fontweight="bold", color=PURPLE)
    if show_pc:
        ax.axvline(d["theta1"], color=DARK, ls="--", lw=1.4); ax.text(d["theta1"], 3, f" PC1: {d['theta1']:.0f}° ({best / total:.0%})", rotation=90, va="bottom", ha="right", fontsize=9)
    ax.set_ylim(0, 105); ax.set_xlim(0, 180); ax.set_xticks([0, 45, 90, 135, 180])
    ax.set_xlabel("angle of the axis (degrees)"); ax.set_ylabel("% of the total variance kept by the projection")
    ax.set_title("② How much variance does each axis keep?")
    # ③ the data seen from the rotated axes
    ax = axes[2]
    ax.scatter(s_u, s_v, s=7, alpha=0.4, color=NS_COLOR, linewidths=0)
    ax.scatter(s_u[d["is_ifn"]], s_v[d["is_ifn"]], s=9, alpha=0.7, color=UP_COLOR, linewidths=0)
    ax.axhline(0, color=PURPLE, lw=1.2); ax.axvline(0, color="#ccc", lw=0.6)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal")
    ax.set_xlabel(f"coordinate on your axis  (variance {np.var(s_u, ddof=1):.2f})"); ax.set_ylabel(f"coordinate on the perpendicular axis  (variance {np.var(s_v, ddof=1):.2f})")
    ax.set_title(f"③ The same samples after rotating by {angle}°\ncorrelation between the new coordinates: {np.corrcoef(s_u, s_v)[0, 1]:+.2f}")
    fig.tight_layout()
    show_fig(fig)
    if abs(((angle - d["theta1"] + 90) % 180) - 90) <= 2:
        msg(f"🎯 Your axis is (almost) <b>PC1</b>: it keeps the maximum variance ({best / total:.0%}) and, in ③, the two new coordinates are no longer correlated: "
            f"<b>PCA is just this rotation</b>. Dropping the perpendicular axis (PC2) would lose only {1 - best / total:.0%} of the variance.", "info")


def exercise_pca_geometry():
    c = dict(pair=dropdown(PAIR_CHOICES, "corr", "520px"), angle=widgets.IntSlider(min=0, max=180, value=0, step=1, continuous_update=False,
                                                                                     layout=_lay("430px"), style=STYLE, readout_format="d"),
             show_proj=checkbox("show the projections (dotted lines)", True, "270px"), show_pc=checkbox("show PC1 and PC2", False, "200px"))
    snap = widgets.Button(description="🎯 jump to PC1", layout=_lay("140px"))
    snap.on_click(lambda _b: setattr(c["angle"], "value", int(round(pair_data(c["pair"].value)["theta1"])) % 181))
    rows = [[field("Choose the pair of genes", c["pair"])],
            [field("Rotate your axis: angle (degrees)", c["angle"], "0° = the x axis (first gene), 90° = the y axis (second gene)"),
             widgets.VBox([widgets.HTML('<div class="cx-label">&nbsp;</div>'), snap]), widgets.VBox([widgets.HTML('<div class="cx-label">&nbsp;</div>'), c["show_proj"]]),
             widgets.VBox([widgets.HTML('<div class="cx-label">&nbsp;</div>'), c["show_pc"]])]]
    ex = Exercise(c, ex8_render, rows, title="Rotate an axis and watch the projection")
    ex.show()
    return ex
