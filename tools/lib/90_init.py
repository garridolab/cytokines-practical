# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  90_init.py: everything slow happens here, ONCE (load the data, normalise, PCA, load the DESeq2 results)
# ══════════════════════════════════════════════════════════════════════════════════════════════════
GENESETS_WANTED = {
    "Monocyte markers": ["LYZ", "CD14", "S100A9", "FCN1"],
    "T-cell markers": ["CD3E", "CD3D"],
    "CD8 T-cell markers": ["CD8A", "CD8B"],
    "NK-cell markers": ["NKG7", "GNLY", "KLRF1"],
    "IFN-γ / STAT1 response": ["STAT1", "IRF1", "GBP1", "IDO1", "CXCL9", "CXCL10"],
    "Type I interferon response": ["ISG15", "IFIT1", "IFIT3", "MX1", "OAS1"],
    "IL-4 / STAT6 response": ["ALOX15", "CCL17", "FCER2"],
    "IL-1β / NF-κB response": ["IL6", "CXCL8", "CCL2"],
    "IL-2 / STAT5 response": ["IL2RA", "GZMB", "CISH"],
    "Housekeeping genes": ["ACTB", "GAPDH", "B2M", "TMSB4X", "UBC"],
}

adata = None          # the AnnData object of the practical (X = raw counts; layers: cpm, log2cpm)
EXPRESSED = None      # True for the genes that pass the expression filter (≥ 10 counts in ≥ 12 samples)
PCA_ALL, PCA_CT, DE = None, {}, None


def init_practical(h5ad_path=None, de_cache_path=None, verbose=True, with_de=True):
    """Load the data and pre-compute everything the exercises need. Takes ~30 s once the files are on disk."""
    global adata, EXPRESSED, PCA_ALL, PCA_CT, DE, GENE_NAMES, GENE_LOOKUP, GENESETS
    t0 = time.time()
    say = (lambda s: print(s, flush=True)) if verbose else (lambda s: None)
    path = get_h5ad(h5ad_path, verbose)
    say("📖 Reading and tidying the data …")
    adata = load_adata(path)
    EXPRESSED = is_expressed(adata.X)
    GENE_NAMES = list(adata.var_names)
    GENE_LOOKUP = {}
    for g in GENE_NAMES:
        GENE_LOOKUP.setdefault(g.upper(), g)
    GENESETS = {k: [g for g in v if g in GENE_LOOKUP.values()] for k, v in GENESETS_WANTED.items()}
    say("🧮 Running the PCAs (all samples, and one per cell type) …")
    PCA_ALL = make_pca(adata, expressed=EXPRESSED)
    PCA_CT = {ct: make_pca(adata, mask=(adata.obs["cell_type"] == ct).values, expressed=EXPRESSED) for ct in CELL_TYPES}
    if with_de:
        say("📦 Loading the differential-expression results …")
        DE = DEResults(adata, de_cache_path or get_de_cache(verbose), verbose)
    say(f"✅ Ready in {time.time() - t0:.0f} s: {adata.n_obs:,} samples × {adata.n_vars:,} genes. Scroll down and start with Exercise 1.")
    return adata
