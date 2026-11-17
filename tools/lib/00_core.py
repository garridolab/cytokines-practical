# ══════════════════════════════════════════════════════════════════════════════════════════════════
#  00_core.py: imports, constants, colours and plot style
#  (The files of tools/lib are concatenated, in alphabetical order, into the hidden setup cell of the notebook.)
# ══════════════════════════════════════════════════════════════════════════════════════════════════
import os, re, io, sys, json, time, math, html, hashlib, difflib, warnings, itertools
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
import ipywidgets as widgets
from IPython.display import display, HTML, clear_output

warnings.filterwarnings("ignore")


# ── the dataset ────────────────────────────────────────────────────────────────────────────────────
GDRIVE_FILE_ID = "1thFoQnNVVqgKvspg5b5k2fBgjk4bEqAn"          # pseudobulk .h5ad on Google Drive ("anyone with the link")
DATA_FILENAME = "pseudobulks_4ct_proteincoding.h5ad"
DE_CACHE_FILENAME = "de_results.npz"                           # pre-computed DESeq2 results (all cytokines x all cell types)
DE_CACHE_URL = "https://raw.githubusercontent.com/garridolab/cytokines-practical/main/cache/" + DE_CACHE_FILENAME

DISPLAY = {}      # display names of a few words that are also data keys (filled in by the Spanish build)


def disp(word):
    return DISPLAY.get(word, word)


CELL_TYPES = ["CD4_T_cell", "CD8_T_cell", "NK_CD56low", "CD14_Mono"]
CT_LABEL = {"CD4_T_cell": "CD4 T cell", "CD8_T_cell": "CD8 T cell", "NK_CD56low": "NK CD56low", "CD14_Mono": "CD14 monocyte"}
CT_SHORT = {"CD4_T_cell": "CD4 T", "CD8_T_cell": "CD8 T", "NK_CD56low": "NK", "CD14_Mono": "Mono"}
CT_COLORS = {"CD4_T_cell": "#2a78d6", "CD8_T_cell": "#008300", "NK_CD56low": "#eda100", "CD14_Mono": "#e87ba4"}
CT_MARKERS = {"CD4_T_cell": "o", "CD8_T_cell": "s", "NK_CD56low": "^", "CD14_Mono": "D"}
DONORS = [f"Donor{i}" for i in range(1, 13)]

# the expression filter of the practical: at least 10 counts in at least 12 samples (12 = the size of one group)
MIN_COUNT, MIN_SAMPLES = 10, 12
N_HVG = 2000                                                    # highly variable genes used for every PCA
PADJ_DEFAULT, LFC_DEFAULT = 0.05, 0.5                           # default thresholds of the volcano plots

# ── colours: "direction" is blue (down) < grey < red (up) ─────────────────────────────────────────────
UP_COLOR, DOWN_COLOR, NS_COLOR, PBS_COLOR = "#e34948", "#2a78d6", "#c3c2b7", "#8a8984"
BLUE, ORANGE, GREEN, YELLOW, PINK, PURPLE = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"
DARK = "#2b2a28"
DIVERGING = LinearSegmentedColormap.from_list("div", [DOWN_COLOR, "#f0efec", UP_COLOR])
SEQUENTIAL = LinearSegmentedColormap.from_list("seq", ["#f4f8fd", "#86b6ef", "#2a78d6", "#0d366b"])
SEQ_POINTS = LinearSegmentedColormap.from_list("seqp", ["#cde2fb", "#2a78d6", "#0d366b"])
DONOR_HUES = [BLUE, ORANGE, GREEN, YELLOW, PINK, PURPLE]       # 12 donors = 6 hues x 2 marker shapes
CYT_HUES = [BLUE, ORANGE, GREEN, YELLOW, PINK, PURPLE, "#008300", "#8a5a00"]
CYT_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "h"]

plt.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 100, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titleweight": "bold", "axes.titlesize": 11, "axes.labelsize": 10, "legend.frameon": False,
    "legend.fontsize": 8.5, "xtick.labelsize": 9, "ytick.labelsize": 9, "figure.autolayout": False,
})
