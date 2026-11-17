# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  10_data.py: download, tidy-up and pre-computation (everything slow happens here, once)
# ══════════════════════════════════════════════════════════════════════════════════════════════════

# ── files ───────────────────────────────────────────────────────────────────────────────────────────
def _search_dirs():
    here = os.getcwd()
    dirs = [here, os.path.join(here, "cache"), os.path.dirname(here), "/content", os.path.expanduser("~")]
    return [d for d in dict.fromkeys(dirs) if os.path.isdir(d)]


def find_file(filename):
    """Path of `filename` if it already exists next to the notebook (or in /content on Colab), else None."""
    for d in _search_dirs():
        p = os.path.join(d, filename)
        if os.path.exists(p):
            return p
    return None


def download_dir():
    return "/content" if os.path.isdir("/content") else os.getcwd()


def get_h5ad(path=None, verbose=True):
    """Path to the pseudobulk .h5ad; downloads it from Google Drive the first time."""
    path = path or find_file(DATA_FILENAME)
    if path and os.path.exists(path):
        return path
    import gdown
    path = os.path.join(download_dir(), DATA_FILENAME)
    if verbose:
        print("⬇️  Downloading the data from Google Drive (≈ 500 MB, about a minute) …")
    out = gdown.download(f"https://drive.google.com/uc?id={GDRIVE_FILE_ID}", path, quiet=not verbose)
    if out is None or not os.path.exists(path):
        raise RuntimeError("The download failed (Google Drive quota?). Wait a minute and run this cell again, or ask the instructor for the file.")
    return path


def get_de_cache(verbose=True):
    """Path to the pre-computed DESeq2 results, or None if they cannot be found/downloaded (the notebook then computes on demand)."""
    p = find_file(DE_CACHE_FILENAME)
    if p:
        return p
    try:
        import urllib.request
        target = os.path.join(download_dir(), DE_CACHE_FILENAME)
        if verbose:
            print("⬇️  Downloading the pre-computed differential-expression results …")
        urllib.request.urlretrieve(DE_CACHE_URL, target)
        return target
    except Exception as err:                                      # offline, or the file is not published (yet)
        if verbose:
            print(f"ℹ️  Pre-computed DE results not available ({type(err).__name__}): they will be computed when needed (≈ 15 s each).")
        return None


# ── tidy-up of the AnnData object ─────────────────────────────────────────────────────────────────
def is_expressed(counts, min_count=MIN_COUNT, min_samples=MIN_SAMPLES):
    """For a samples × genes matrix of raw counts: True for the genes with at least `min_count` counts in at least `min_samples` samples."""
    return (np.asarray(counts) >= min_count).sum(axis=0) >= min_samples


def load_adata(path):
    """Read the .h5ad and return the tidy AnnData of the practical:
    X = raw counts (dense), obs = sample metadata, var = gene metadata, layers: `cpm` and `log2cpm` (normalised versions of X)."""
    raw = ad.read_h5ad(path)
    # 1) the PBS control was measured in 6 wells and the file stores each well *plus* their total: keep the total only
    is_pbs_well = ((raw.obs["cytokine"] == "PBS") & raw.obs["bc1_well"].notna()).values
    raw = raw[~is_pbs_well]
    # 2) sample metadata, with the order of the categories fixed, and a unique readable name:  donor|cell type|cytokine
    obs = raw.obs[["donor", "cytokine", "treatment", "cell_type", "psbulk_n_cells"]].copy()
    cytokines = ["PBS"] + sorted(c for c in obs["cytokine"].astype(str).unique() if c != "PBS")
    obs["donor"] = pd.Categorical(obs["donor"].astype(str), categories=DONORS)
    obs["cytokine"] = pd.Categorical(obs["cytokine"].astype(str), categories=cytokines)
    obs["cell_type"] = pd.Categorical(obs["cell_type"].astype(str), categories=CELL_TYPES)
    obs["treatment"] = pd.Categorical(np.where(obs["cytokine"] == "PBS", "control", "treated"), categories=["control", "treated"])
    obs["psbulk_n_cells"] = obs["psbulk_n_cells"].astype(int)
    obs.index = (obs["donor"].astype(str) + "|" + obs["cell_type"].astype(str) + "|" + obs["cytokine"].astype(str)).values
    assert obs.index.is_unique
    # 3) the matrix: raw counts, as a plain dense table of numbers
    X = raw.X.toarray() if hasattr(raw.X, "toarray") else np.asarray(raw.X)
    X = np.ascontiguousarray(X, dtype=np.float32)
    var = pd.DataFrame(index=pd.Index(raw.var_names.astype(str), name=None))
    del raw
    # 4) our own columns: library size (obs) and total counts / detection (var)
    obs["library_size"] = X.sum(axis=1, dtype=np.float64).astype(np.int64)
    var["total_counts"] = X.sum(axis=0, dtype=np.float64).astype(np.int64)
    var["n_samples_detected"] = (X > 0).sum(axis=0).astype(int)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    # 5) normalised versions of the matrix, kept in `layers` (more matrices with the same shape as X)
    lib = obs["library_size"].values.astype(np.float64)[:, None]
    cpm = (X / lib * 1e6).astype(np.float32)
    adata.layers["cpm"] = cpm
    adata.layers["log2cpm"] = np.log2(cpm + 1).astype(np.float32)
    return adata


# ── PCA ────────────────────────────────────────────────────────────────────────────────────────────
def pca_svd(M, n_components=10):
    """PCA of a samples × genes matrix: centre → SVD → project. Returns scores, % of variance explained, loadings (genes × PCs)."""
    Mc = M.astype(np.float64) - M.mean(axis=0)                    # 1) centre every gene
    U, S, Vt = np.linalg.svd(Mc, full_matrices=False)             # 2) SVD
    for k in range(min(n_components, len(S))):                    # the sign of a PC is arbitrary: make it reproducible
        j = np.argmax(np.abs(Vt[k]))
        if Vt[k, j] < 0:
            Vt[k], U[:, k] = -Vt[k], -U[:, k]
    scores = (U * S)[:, :n_components]                            # 3) projection: coordinates of the samples on the PCs
    var_pct = (S ** 2 / np.sum(S ** 2) * 100)[:n_components]
    return scores, var_pct, Vt[:n_components].T


def r2_categorical(y, groups):
    """Fraction of the variance of y explained by group membership (one-way ANOVA R²)."""
    y, g = np.asarray(y, float), np.asarray(groups).astype(str)
    ss_tot = ((y - y.mean()) ** 2).sum()
    codes, inv = np.unique(g, return_inverse=True)
    sums, counts = np.bincount(inv, weights=y), np.bincount(inv)
    return float((counts * (sums / counts - y.mean()) ** 2).sum() / ss_tot)


def pc_r2_table(scores, obs, n_pc=10):
    """R² of every PC with every sample variable (categories: ANOVA R²; numeric: squared Pearson correlation)."""
    loglib = np.log10(obs["library_size"].astype(float).values)
    logn = np.log10(obs["psbulk_n_cells"].astype(float).values)
    cols = {
        "cell type": lambda y: r2_categorical(y, obs["cell_type"]),
        "donor": lambda y: r2_categorical(y, obs["donor"]),
        "cytokine": lambda y: r2_categorical(y, obs["cytokine"]),
        "control vs treated": lambda y: r2_categorical(y, obs["treatment"]),
        "library size": lambda y: np.corrcoef(y, loglib)[0, 1] ** 2,
        "number of cells": lambda y: np.corrcoef(y, logn)[0, 1] ** 2,
    }
    n_pc = min(n_pc, scores.shape[1])
    return pd.DataFrame({k: [f(scores[:, i]) if np.std(scores[:, i]) > 0 else 0.0 for i in range(n_pc)] for k, f in cols.items()},
                        index=[f"PC{i + 1}" for i in range(n_pc)])


def make_pca(adata, mask=None, n_hvg=N_HVG, expressed=None):
    """PCA on the `n_hvg` most variable genes (log2 CPM + 1) of the selected samples. Returns a dict with everything the plots need."""
    mask = np.ones(adata.n_obs, bool) if mask is None else np.asarray(mask)
    L = adata.layers["log2cpm"][mask]
    ok = np.ones(adata.n_vars, bool) if expressed is None else np.asarray(expressed)
    cand = np.flatnonzero(ok)
    var_g = L[:, cand].var(axis=0)
    hvg = cand[np.argsort(var_g)[::-1][:n_hvg]]
    scores, var_pct, loadings = pca_svd(L[:, hvg])
    obs = adata.obs[mask]
    return dict(mask=mask, obs=obs, scores=scores, var=var_pct, loadings=loadings, hvg=hvg,
                hvg_names=adata.var_names[hvg].to_numpy(), r2=pc_r2_table(scores, obs))


def ensure_pydeseq2():
    """pydeseq2 is only needed if the pre-computed results are missing: install it the first time it is needed."""
    import importlib.util, subprocess
    if importlib.util.find_spec("pydeseq2") is None:
        print("📦 Installing pydeseq2 …")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "pydeseq2>=0.5.4"])
        importlib.invalidate_caches()


# ── differential expression (DESeq2) ──────────────────────────────────────────────────────────────────
def de_gene_universe(adata, expressed=None):
    """The genes that take part in the DE exercises: the immune list ∩ the genes expressed in the whole data set (≈ 500)."""
    expressed = is_expressed(adata.X) if expressed is None else expressed
    return [g for g in adata.var_names[expressed] if g in set(IMMUNE_GENES)]


def run_deseq2(adata, cell_type, treated, reference="PBS", n_cpus=1, genes=None):
    """DESeq2 (pydeseq2): `treated` vs `reference` in one cell type, 12 donors, design `~ donor + cytokine`, restricted to `genes`
    (default: the immune genes, see de_gene_universe). One row per gene of `adata` (NaN for the genes that were not tested).
    This is what produced the cached results."""
    genes = de_gene_universe(adata) if genes is None else genes
    ensure_pydeseq2()
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    names = [f"{d}|{cell_type}|{c}" for c in (reference, treated) for d in DONORS]          # 12 reference + 12 treated samples
    rows = adata.obs_names.get_indexer(names)
    counts = pd.DataFrame(np.rint(adata.X[rows][:, adata.var_names.get_indexer(genes)]).astype(int), index=names, columns=genes)   # RAW counts, integers
    counts = counts.loc[:, is_expressed(counts.values)]                                      # expression filter on these 24 samples
    meta = pd.DataFrame({"donor": DONORS * 2, "cytokine": [reference] * 12 + [treated] * 12}, index=names)
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~ donor + cytokine", quiet=True, n_cpus=n_cpus)
    dds.deseq2()                                                                             # size factors → dispersions → model fit
    stat = DeseqStats(dds, contrast=["cytokine", treated, reference], quiet=True, n_cpus=n_cpus)
    stat.summary()                                                                           # Wald test + Benjamini-Hochberg
    res = stat.results_df.rename(columns={"log2FoldChange": "log2fc", "lfcSE": "lfc_se", "pvalue": "pval", "baseMean": "base_mean"})
    return res[["base_mean", "log2fc", "lfc_se", "pval", "padj"]].reindex(adata.var_names)


CACHE_FIELDS = ("lfc", "nlp_adj")


def de_to_arrays(res):
    """Compact (float16) version of a DESeq2 table, for the cache: log2 fold change and -log10(adjusted p-value)."""
    f16 = lambda v: np.asarray(v, dtype=np.float64).astype(np.float16)
    return dict(lfc=f16(res["log2fc"].clip(-60, 60)), nlp_adj=f16(-np.log10(res["padj"].clip(lower=1e-300))))


class DEResults:
    """All the differential-expression results (cytokine vs PBS, per cell type), pre-computed once with DESeq2.
    `get(cell_type, cytokine)` returns a table with one row per gene; if the cache file is missing, it is computed on demand."""

    def __init__(self, adata, cache_path=None, verbose=True):
        self.adata, self.live, self._memo = adata, {}, {}
        self.arrays = None
        if cache_path and os.path.exists(cache_path):
            z = np.load(cache_path, allow_pickle=False)
            genes = pd.Index(z["genes"].astype(str))
            take = genes.get_indexer(adata.var_names)                              # re-order the genes like `adata`
            self.cell_types = [str(c) for c in z["cell_types"]]
            self.cytokines = [str(c) for c in z["cytokines"]]
            self.arrays = {k: self._reindex(z[k], take) for k in CACHE_FIELDS}
            self.meta = json.loads(str(z["meta"])) if "meta" in z.files else {}
            if verbose:
                print(f"✅ Pre-computed DESeq2 results: {len(self.cell_types)} cell types × {len(self.cytokines)} cytokines")
        else:
            self.cell_types, self.cytokines, self.meta = list(CELL_TYPES), [c for c in adata.obs["cytokine"].cat.categories if c != "PBS"], {}

    @staticmethod
    def _reindex(a, take):
        out = np.full(a.shape[:2] + (len(take),), np.nan, dtype=np.float16)
        ok = take >= 0
        out[..., ok] = a[..., take[ok]]
        return out

    @property
    def cached(self):
        return self.arrays is not None

    def get(self, cell_type, cytokine):
        key = (cell_type, cytokine)
        if key in self._memo:
            return self._memo[key]
        if self.cached and cell_type in self.cell_types and cytokine in self.cytokines:
            i, j = self.cell_types.index(cell_type), self.cytokines.index(cytokine)
            nlp_adj = self.arrays["nlp_adj"][i, j].astype(np.float64)
            res = pd.DataFrame({"log2fc": self.arrays["lfc"][i, j].astype(np.float32), "padj": 10.0 ** -nlp_adj,
                                "neglog10_padj": nlp_adj}, index=self.adata.var_names)
        else:                                                                       # not in the cache: run DESeq2 now (≈ 15 s)
            res = run_deseq2(self.adata, cell_type, cytokine, n_cpus=2)[["log2fc", "padj"]].copy()
            res["neglog10_padj"] = -np.log10(res["padj"].clip(lower=1e-300))
        res["tested"] = res["padj"].notna()
        self._memo[key] = res
        return res

    def lfc_matrix(self, cell_type):
        """cytokines × genes matrix of log2 fold changes in one cell type (cache only)."""
        return self.arrays["lfc"][self.cell_types.index(cell_type)].astype(np.float32) if self.cached else None


def classify_de(res, padj_thresh=PADJ_DEFAULT, lfc_thresh=LFC_DEFAULT):
    """'up' / 'down' / 'ns' for every gene of a DE table (`ns` also for the genes that were not tested)."""
    sig = res["tested"] & (res["padj"] < padj_thresh)
    direction = pd.Series("ns", index=res.index)
    direction[sig & (res["log2fc"] > lfc_thresh)] = "up"
    direction[sig & (res["log2fc"] < -lfc_thresh)] = "down"
    return direction
