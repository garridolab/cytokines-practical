#!/usr/bin/env python
"""Pre-computes the differential expression (DESeq2 via pydeseq2) of EVERY cytokine vs PBS in EVERY cell type
(87 cytokines x 4 cell types = 348 comparisons, design `~ donor + cytokine`, paired by donor) and stores it in a compact file,
so that the students' notebook never has to run DESeq2.

Usage:   python tools/precompute_de.py [--workers 16] [--h5ad pseudobulks_4ct_proteincoding.h5ad] [--out cache/de_results.npz]
Needs:   pydeseq2 >= 0.5.4, anndata, the .h5ad file.  Takes ≈ 5-10 minutes with 16 workers (≈ 11 s per comparison and core).
Output:  cache/de_results.npz   (float16 arrays [cell type, cytokine, gene]: lfc, nlp_adj; see `de_to_arrays` in tools/lib/10_data.py)
"""
import argparse, json, os, sys, time, datetime, multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from libsource import load_library

L = load_library()
ADATA = None
GENES = None


def _work(job):
    ct, cyt = job
    t0 = time.time()
    res = L["run_deseq2"](ADATA, ct, cyt, n_cpus=1, genes=GENES)
    return ct, cyt, L["de_to_arrays"](res.loc[GENES]), time.time() - t0


def main():
    global ADATA, GENES
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--h5ad", default=str(ROOT / L["DATA_FILENAME"]))
    ap.add_argument("--out", default=str(ROOT / "cache" / L["DE_CACHE_FILENAME"]))
    ap.add_argument("--only", default=None, help="debug: comma-separated cytokines to compute")
    args = ap.parse_args()

    print("loading data …", flush=True)
    ADATA = L["load_adata"](args.h5ad)
    GENES = L["de_gene_universe"](ADATA)                                              # ≈ 500 immune genes
    print(f"{len(GENES)} genes in the DE analysis", flush=True)
    cell_types = list(L["CELL_TYPES"])
    cytokines = [c for c in ADATA.obs["cytokine"].cat.categories if c != "PBS"]
    if args.only:
        cytokines = [c for c in cytokines if c in args.only.split(",")]
    jobs = [(ct, cyt) for ct in cell_types for cyt in cytokines]
    n_genes = len(GENES)
    import numpy as np
    arrays = {k: np.full((len(cell_types), len(cytokines), n_genes), np.nan, dtype=np.float16) for k in L["CACHE_FIELDS"]}
    t0 = time.time()
    ctx = mp.get_context("fork")                                                    # workers inherit the loaded data
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx) as ex:
        futs = [ex.submit(_work, j) for j in jobs]
        for n, f in enumerate(as_completed(futs), 1):
            ct, cyt, arr, dt = f.result()
            i, j = cell_types.index(ct), cytokines.index(cyt)
            for k in arrays:
                arrays[k][i, j] = arr[k]
            if n % 10 == 0 or n == len(jobs):
                el = time.time() - t0
                print(f"[{n:3d}/{len(jobs)}] {ct:11s} {cyt:15s} {dt:4.1f}s   elapsed {el / 60:.1f} min, ETA {(el / n) * (len(jobs) - n) / 60:.1f} min", flush=True)

    import pydeseq2
    meta = dict(created=datetime.date.today().isoformat(), pydeseq2=getattr(pydeseq2, "__version__", "?"), design="~ donor + cytokine",
                reference="PBS", filter=f">= {L['MIN_COUNT']} counts in >= {L['MIN_SAMPLES']} of the 24 samples of the comparison",
                n_comparisons=len(jobs), padj="Benjamini-Hochberg, per comparison", lfc_storage="float16", gene_set=f"{len(GENES)} expressed immune genes (cytokine signalling, T-cell activation, interferon response)")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, genes=np.array(GENES, dtype=str), cell_types=np.array(cell_types), cytokines=np.array(cytokines),
                        meta=np.array(json.dumps(meta)), **arrays)
    print(f"✅ saved {args.out}  ({os.path.getsize(args.out) / 1e6:.1f} MB) in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
