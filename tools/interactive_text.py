"""Markdown narrative of the interactive notebook (one dictionary per language, same keys).

Every number quoted in an answer was checked against the data with the default settings of the controls
(see tools/smoke_test_interactive.py and the 'numbers' section of tools/build_interactive_notebooks.py).
If you change a default or the data, re-check them.
"""

COLAB_URL = "https://colab.research.google.com/github/garridolab/cytokines-practical/blob/main/{notebook}"
BADGE = "[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({url})"

EN = {}

# ──────────────────────────────────────────────────────────────────────────────────────────────────
EN["title"] = """# 🧬 Human Cytokine Dictionary: an interactive transcriptomics practical

{badge}

**Based on:** Oesinghaus *et al.* (2025), *A single-cell cytokine dictionary of human peripheral blood*, bioRxiv
([doi:10.64898/2025.12.12.693897](https://doi.org/10.64898/2025.12.12.693897)).

In the paper, ~9.7 million human blood cells (PBMCs) from **12 donors** were stimulated with **90 cytokines** (or with a PBS control)
and read out by single-cell RNA-seq. Here you work with a **pseudobulk** version of the data: for every *donor × cytokine × cell type*,
the counts of all the cells were **added up** into one bulk-like expression profile. We keep **four cell types** (CD4 T cells, CD8 T cells, NK cells, CD14 monocytes),
**87 cytokines** and the **PBS control**.

### 🎛 You will not write any code
This is an **interactive** practical. Every exercise gives you a few **controls** (menus, sliders, boxes) and a **figure that redraws by itself** each time you change one.
Your job is to **play, observe and interpret**: what changes? why?

| Part | Exercises | Time |
|:--:|---|:--:|
| **A · Python warm-up** | **1** Variables · **2** Basic plots | 12 min |
| **B · Meet the data** | **3** The AnnData object · **4** Distributions with boxplots | 18 min |
| **C · Normalization** | **5** Raw counts vs CPM · **6** The logarithm | 18 min |
| **D · Patterns in the data** | **7** Heatmaps · **8** PCA: the idea · **9** PCA: all samples · **10** PCA: one cell type | 38 min |
| **E · Differential expression** | **11** Volcano plots · **12** Comparing two responses | 22 min |
| | Wrap-up and quiz | 7 min |

### How to use this notebook

| | |
|:--:|---|
| ▶ | **Run a cell** with the ▶ button on its left (or click it and press **Shift + Enter**). Run the **setup cell** below first, and then each exercise cell *once*: the controls appear. |
| 🎛 | **Change a control and wait a second**: the figures and tables update automatically. You never need to run the exercise cell again. |
| 🖱 | In a list, **Ctrl/⌘ + click** selects several items. In text boxes, press **Enter** after typing. |
| ❓ | Each exercise ends with **questions**. Think (or discuss with your neighbour) *before* opening the 💡 answer. |
| 🔁 | If something looks odd, re-run the exercise cell. If the notebook was restarted, run the setup cell again. |

> 💡 In Google Colab the code of the grey cells is hidden on purpose. If you are curious, double-click the cell title to see it.
"""

EN["setup_md"] = """---
# ⚙️ Setup: run this cell first

It installs what is missing, downloads the data (≈ 500 MB, about a minute on Colab), normalizes it, runs the PCAs and loads the differential-expression
results that were computed beforehand with DESeq2. **You only run it once** (and again if the notebook restarts). Wait until it prints **✅ Ready**."""

# ──────────────────────────────────────────────────────────────────────────────────────────────────
EN["partA"] = """---
# Part A · Python warm-up
Two short exercises that need no data: what a **variable** is, and how to **look at numbers** with plots."""

EN["ex1"] = """## Exercise 1 · Variables: giving names to values
⏱ *5 min*

A **variable** is a **name** that stores a **value**, like a labelled tube in the freezer. You create one with `=`: the name on the left, the value on the right.
Python needs to know *what kind* of value it is (its **type**): **text** (`str`, always between quotes), **whole numbers** (`int`), **decimal numbers** (`float`), **True/False** (`bool`) or a **list** of values.
`print(...)` shows a value on the screen, and variables can be combined to create new ones.

**What to do:** change the controls and watch how the Python code, the printed output and the table of variables change."""

EN["ex1_try"] = """#### ❓ Questions
1. Change `message` to your own name. Which lines of the code panel and of the output change?
2. Move the `n_donors` slider. Which printed lines change and which do not? Why?
3. In the table of variables, which type is `library_size`? Why is it not a whole number?

<details>
<summary>💡 Answers</summary>

1. The line `message = '...'` and the first printed line. Text values are written **between quotes**; the variable keeps whatever you type.
2. Only the lines that *use* `n_donors`: `Donors: …` and `… was added to … donors`. A change in a variable only affects the places where its name is used.
3. It is a `float` (decimal): `library_size = n_cells * umis_per_cell` multiplies a whole number by a decimal number, and the result is a decimal. If both were whole numbers (`int`), the result would be an `int`. (`library_size`, the total number of counts of a sample, is a concept we will meet again in Exercise 3.)

</details>"""

EN["ex2"] = """## Exercise 2 · Looking at numbers: basic plots
⏱ *7 min*

A plot turns many numbers into one picture. Different plots answer different questions:

| Plot | Answers the question… |
|---|---|
| 📊 **Histogram** | *How are the values distributed?* (how many donors have each value) |
| 📦 **Box plot** | *Where is the middle and how spread out are the values?* Line = median, box = the middle 50 %, whiskers = the range |
| 📶 **Bar chart** | *What is the average of each group?* (it hides the individual values!) |
| 🔵 **Scatter plot** | *Are two measurements related?* One point per donor |
| 📈 **Line plot** | *How does a value change over time?* |

**What to do:** the data are **simulated**: one gene measured in several donors, without (control) and with (treated) a cytokine. Pick a plot type and play with the numbers."""

EN["ex2_try"] = """#### ❓ Questions
1. In the **histogram**, set the effect of the cytokine to 0 and then to +2. What moves? What happens if you also increase the variability between donors (SD)?
2. Use few donors (5–10) and press **🎲 new random draw** several times (histogram or box plot). Do the two groups always look the same? What does this tell you about experiments with few donors?
3. In the **scatter plot**, use an effect of +1 and compare a similarity of 0.1 with one of 0.95. When is the effect easier to see? Why would it help to compare every donor with *itself*?

<details>
<summary>💡 Answers</summary>

1. With effect 0 the two histograms overlap completely; with +2 the treated one shifts to the right by 2 units (2 units of log₂ = **4 times** more expression). A larger SD makes both histograms wider, so they overlap more even if the effect is the same: the effect is harder to see when donors differ a lot.
2. No: with few donors, chance alone changes the means and the shapes, and you can see "differences" (or miss real ones) by bad luck. More donors give stable pictures. This is why experiments need replicates.
3. With similarity 0.95 all the points sit in a narrow band **above** the dashed line (`y = x`): *every* donor goes up by about the same amount, even though donors start from very different levels. With 0.1 the cloud is round and the shift is hidden by the differences between donors. Comparing each donor with itself (a **paired** design, Exercises 11–12) removes the donor-to-donor differences.

</details>"""

# ──────────────────────────────────────────────────────────────────────────────────────────────────
EN["partB"] = """---
# Part B · Meet the data"""

EN["ex3"] = """## Exercise 3 · The AnnData object: samples, genes and their metadata
⏱ *10 min*

Our data live in an **AnnData** object called `adata`, which keeps several tables **aligned**:

| Component | What it holds | Shape |
|---|---|---|
| `X` | the **numbers**: raw counts (how many mRNA molecules of each gene were detected) | samples × genes |
| `obs` | information about the **samples** (rows of `X`): donor, cytokine, cell type… | one row per sample |
| `var` | information about the **genes** (columns of `X`) | one row per gene |
| `layers` | other versions of `X` with the same shape (here, *normalized* versions we will build in Exercises 5–6) | samples × genes |

The complete object has **4,224 samples × 18,120 genes**: 12 donors × 4 cell types × 88 conditions (87 cytokines + PBS).
A **sample** is a pseudobulk: all the cells of one donor, one cell type and one condition, added together. Its name is `donor|cell type|cytokine`.
Selecting works like a table: `adata[samples, genes]`, **samples first, genes second**.

**What to do:** choose samples and genes (the defaults are the 4 PBS controls of Donor 1, plus its monocytes treated with IFN-γ). The **subset** you build is shown in the three tabs, and its place in the full matrix on the right."""

EN["ex3_try"] = """#### ❓ Questions
1. In the default subset, which sample has the highest count of `LYZ`, of `CD3E` and of `NKG7`? Does it match what you know about these genes?
2. Open the **obs** tab: what do `psbulk_n_cells` and `library_size` tell you? Compare the two monocyte samples (PBS and IFN-γ).
3. `STAT1` is switched on by IFN-γ. Compare its count in the two monocyte samples. Is it higher with IFN-γ? Is that what you expected?

<details>
<summary>💡 Answers</summary>

1. `LYZ` is highest in the **monocytes** (911,945 counts), `CD3E` in the **CD4 T cells** (20,361; CD8 T cells are second with 3,100) and `NKG7` in the **NK cells** (2,392; then CD8 T cells, which are also cytotoxic): they are cell-type markers.
2. `psbulk_n_cells` is the **number of cells that were added up** into the sample, and `library_size` the **total number of counts** (of all genes). The PBS monocyte sample contains 17,373 cells and 147.9 million counts; the IFN-γ one only 2,556 cells and 21.4 million counts: about **7 times less material**. (The PBS control is the sum of 6 wells.)
3. **No: it is ~9 times *lower* with IFN-γ** (73,853 vs 639,289), although IFN-γ induces `STAT1`. The raw count mixes **biology** (is the gene switched on?) with the **amount of material** (how many cells, how deep). That is the problem we solve in Part C.

</details>"""

EN["ex4"] = """## Exercise 4 · Distributions of raw counts: boxplots
⏱ *8 min*

A **box plot** summarizes many values: the **line** is the median, the **box** contains the middle 50 % of the values and the **whiskers** show the range. Each **dot** is one **pseudobulk sample**
(a donor with one cytokine in one cell type); the *box* is therefore the distribution of the gene across all samples of a group.
The right-hand panel shows the **library size** (total counts) of the same samples.

**What to do:** the default shows four genes with different roles in the four cell types (all 88 conditions, 1,056 samples per cell type). Change the genes, the groups or the cell types.
Try "control (PBS) vs treated" to compare the two kinds of samples."""

EN["ex4_try"] = """#### ❓ Questions
1. In the default plot, how many orders of magnitude separate the most and the least abundant gene in monocytes? Do `LYZ`, `CD3E` and `NKG7` behave as cell-type markers?
2. `ACTB` is a housekeeping gene (always on), yet its counts differ between cell types. Compare with the right panel: is there a link?
3. Group the boxes by **control (PBS) vs treated**. What do you see for *every* gene? Why?

<details>
<summary>💡 Answers</summary>

1. In monocytes the median raw count goes from 4 (`NKG7`) to 25,725 (`STAT1`): **almost 4 orders of magnitude**. `LYZ` is highest in monocytes (median 30,058 vs ≤ 297 in the other cells), `CD3E` in CD4 and CD8 T cells (1,262 and 663) and `NKG7` in NK cells (134; then CD8 T cells, 65).
2. Yes: the **library size** differs between cell types (medians: NK 1.2 M, CD8 T 3.8 M, monocytes 7.6 M, CD4 T 8.7 M), and `ACTB` counts follow it largely (NK lowest, CD8 T next). Monocytes are bigger cells with more mRNA, so they have even more `ACTB` per count. Counts depend on how much material each cell type contributes.
3. The **PBS controls are higher for every gene**, including `ACTB`: their library size is ~5 times larger (median 23 M vs 4.4 M) because each PBS sample is the sum of 6 wells. This difference has nothing to do with the biology of each gene, it is the amount of material, which we correct in the next exercise.

</details>"""

# ──────────────────────────────────────────────────────────────────────────────────────────────────
EN["partC"] = """---
# Part C · Normalization: making samples comparable"""

EN["ex5"] = """## Exercise 5 · Raw counts vs CPM: why do we normalize?
⏱ *10 min*

A raw count depends on **two things**: how much the gene is expressed (biology) and how much material we have (how many cells, how deeply sequenced; technical).
The **library size** is the total number of counts of a sample. To remove the effect of depth we express every gene as a fraction of its sample:

**CPM** (counts per million) = count of the gene ÷ library size of the sample × 1,000,000

**What to do:**
* **Toy example** (simulated): the same donor with a deep sample (PBS) and a shallow one (with a cytokine). Only gene B really changes. Move the sliders.
* **Real data, PBS vs cytokine**: one line per donor. The default is the housekeeping gene `ACTB` in monocytes with IFN-γ. Change the gene, the cell type or the cytokine.
* **Real data, all samples**: does the count of the gene follow the library size?"""

EN["ex5_try"] = """#### ❓ Questions
1. **Toy example**: with a depth of 20 % and a true fold change of ×4, what do the raw counts say about gene B? And the CPM? What if the depth is 100 %?
2. **Real data**, `ACTB` in monocytes (PBS vs IFN-γ): how do the library size, the raw `ACTB` counts and the `ACTB` CPM change? What do you expect from a housekeeping gene?
3. Change the gene to `IDO1` and then to `STAT1` (both induced by IFN-γ). What do raw counts and CPM say? Then look at the **all samples** view: how does the correlation with the library size change from raw to CPM?

<details>
<summary>💡 Answers</summary>

1. **Raw counts: gene B looks slightly *down* (×0.78)** and every other gene looks ~5 times lower, simply because sample 2 is 5 times shallower. **CPM: genes A, C and D are unchanged (×0.97) and gene B is ×3.9**, close to the true ×4. With a depth of 100 % both representations agree: normalization only matters when depths differ.
2. The library size is **6 times lower** with IFN-γ (×0.16) and the raw `ACTB` count follows it (×0.16); after CPM `ACTB` is **unchanged (×0.95)**, as expected for a housekeeping gene. In the all-samples view the correlation of `ACTB` with the library size is **r = +0.84 in raw counts and −0.25 in CPM**: the depth effect is gone.
3. `IDO1`: raw counts say **no change (×1.1)**, CPM says **×6.7 up**. `STAT1`: raw counts say **×4.5 *down***, CPM says **×1.35 up**. Without normalization we would miss the induction of `IDO1` and conclude the opposite for `STAT1`.

</details>"""

EN["ex6"] = """## Exercise 6 · The logarithm: taming the scale
⏱ *8 min*

CPM fixes the differences between **samples**, but not between **genes**: some genes have <1 CPM and others >10,000, the values are very **skewed**, and the **variance grows with the mean**.
Many analyses (PCA, tests…) are driven by variance, so a handful of abundant genes would decide everything. The **logarithm** compresses large values and turns ratios into differences
(with log₂: doubling = +1, halving = −1). Because log(0) does not exist, we first add a small **pseudocount** (usually 1): **log₂(CPM + 1)**.

**What to do:** the top panels show one gene in one cell type; the bottom panels show **all genes** (variance against mean). Change the **base** of the logarithm and the **pseudocount**."""

EN["ex6_try"] = """#### ❓ Questions
1. Compare ① and ② for `CXCL10` in monocytes. What happens to the shape of the distribution and to its variance?
2. Compare ③ and ④: which share of the total variance do the 20 most variable genes hold, before and after the logarithm? What would that mean for a PCA?
3. Change the base to 10 and to *e*. What changes and what does not? And with the pseudocount at 0 and at 10?

<details>
<summary>💡 Answers</summary>

1. In CPM the distribution is **skewed** (long tail: mean 1,297, median 634, skewness 1.6) and its variance is huge (2,350,366). After log₂ it is much more symmetric (skewness −0.3) and the variance is 4.7.
2. In CPM **20 genes hold 41 % of the whole variance**; after log₂(CPM + 1) only **3 %**. Without the logarithm, a PCA would be decided by a few abundant genes (and would mostly "see" them); after it, thousands of genes contribute.
3. The **shape does not change**, only the numbers on the axis (log₂ x = 3.32 × log₁₀ x): the variances are multiplied by 0.48 (ln) or 0.09 (log₁₀). The base is a matter of convention; we use 2 because a difference of 1 means "twice as much". With pseudocount **0**, zeros cannot be transformed and are dropped. With **10**, low values are squashed together. 1 is the usual compromise.

</details>"""

# ──────────────────────────────────────────────────────────────────────────────────────────────────
EN["partD"] = """---
# Part D · Finding patterns in the data
From now on we use the **normalized** values: log₂(CPM + 1)."""

EN["ex7"] = """## Exercise 7 · Heatmaps: genes × samples
⏱ *10 min*

A **heatmap** shows many genes (rows) in many samples (columns) at once. To compare patterns, every gene (row) is **scaled across the samples** (z-score): **red** = above the average of *that* gene, **blue** = below.
The two bars on top show the **cell type** and the **library size** of each column (dark = deep).

**What to do:** the default shows cell-type markers and cytokine-response genes in the four cell types with PBS and four cytokines (12 donors each), with the **raw counts** on the left and the **normalized** values on the right.
Switch between the views, choose other **gene groups**, add your own genes (e.g. `TNF`, `IL1B`, `IL13RA1`) and select other **cytokines** (Ctrl/⌘ + click)."""

EN["ex7_try"] = """#### ❓ Questions
1. Compare BEFORE and AFTER. Which columns are the brightest in the raw heatmap, for almost every gene? Can you read the cell types in each panel?
2. In the normalized heatmap, which cytokine switches on which block of genes, and in which cell type?
3. Select **IFN-beta** and the group *Type I interferon response*. In which cell types do these genes respond? Then add `IL13RA1` and `IL4R` as extra genes: which cell type expresses `IL13RA1`? (keep this in mind for Exercise 12)

<details>
<summary>💡 Answers</summary>

1. **Before**, the brightest columns are the **PBS** ones (dark in the library-size bar), whatever the gene: the picture shows *how many counts each sample has*, and the cell types are hard to read. **After**, the PBS glow is gone and every cell type shows its own markers: `LYZ`, `CD14`, `S100A9`, `FCN1` in monocytes; `CD3E`, `CD3D` in T cells; `CD8A`, `CD8B` in CD8 T cells; `NKG7`, `GNLY`, `KLRF1` in NK cells.
2. **IL-4** → `ALOX15`, `CCL17`, `FCER2` almost only in **monocytes**. **IFN-γ** → `GBP1`, `IDO1`, `CXCL9`, `CXCL10`, strongest in **monocytes**. **IL-1β** → `IL6`, `CXCL8`, `CCL2` in **monocytes**. **IL-2** → `IL2RA`, `GZMB`, `CISH` in **NK and T cells**, the lymphocytes that carry the IL-2 receptor.
3. The interferon-stimulated genes (`ISG15`, `IFIT1`, `MX1`…) go up **in all four cell types** with IFN-β (log₂ fold change between 2 and 3), because the type I interferon receptor is present in all of them. `IL13RA1` (a subunit of the IL-13 receptor) is expressed in **monocytes** (~256 CPM) and practically absent from T and NK cells (~2 CPM).

</details>

> ⚠️ **z-scores are for visualization.** They hide the absolute level of each gene and exaggerate genes that barely change. Never compute fold changes from them."""

EN["ex8"] = """## Exercise 8 · PCA: the idea, with just two genes
⏱ *10 min*

**PCA** (principal component analysis) summarizes many genes in a few new axes. To see how, take just **two genes**: each sample is a point in a plane (*x* = gene 1, *y* = gene 2).
1. Draw an **axis** (a line through the centre of the cloud) and **project** every sample onto it (drop a perpendicular from the point to the line).
2. The **spread (variance)** of the projected points tells how much information that axis keeps.
3. **Rotate** the axis until the variance is **maximal**: that axis is **PC1**. The perpendicular one is **PC2**.

PCA is *just a rotation* of the axes (no information is lost); **dimensionality reduction** appears when we keep only the first PCs. With thousands of genes the recipe is exactly the same.

**What to do:** there are two fixed examples in monocytes: a **correlated** pair of genes (`ISG15` and `MX1`, both switched on by interferon) and an **uncorrelated** pair (`IDO1` and `IL1B`, two different pathways).
Rotate the axis with the slider and watch panels ①②③. Tick *show PC1 and PC2*, or press *jump to PC1*."""

EN["ex8_try"] = """#### ❓ Questions
1. **Correlated pair**: at which angle does the axis keep the most variance? How much? What does the axis look like compared with the cloud of points?
2. **Uncorrelated pair**: how much does the kept variance change with the angle? Could one axis summarize both genes?
3. In panel ③, at the angle of PC1, what is the correlation between the two new coordinates? If we drop PC2, how much variance do we lose, in each pair?

<details>
<summary>💡 Answers</summary>

1. At about **36°** the axis keeps **98 %** of the total variance: it runs along the long direction of the elongated cloud (the two genes rise and fall together, r = 0.95). One number per sample (the position along PC1) can replace the two genes.
2. Very little: between **45 % and 55 %** whatever the angle (the best axis, PC1, keeps 55 %). The cloud has no preferred direction, so no single axis summarizes both genes: we would lose almost half of the variance.
3. **About 0**: PC1 and PC2 are uncorrelated by construction. Dropping PC2 loses only **2 %** of the variance for the correlated pair but **45 %** for the uncorrelated pair. PCA compresses well *when genes are correlated*: with thousands of genes that respond together (all the interferon genes, all the monocyte genes…) a few PCs summarize most of the variation.

</details>"""

EN["ex9"] = """## Exercise 9 · PCA of all the samples: what drives the variation?
⏱ *8 min*

Each of the 4,224 samples is described by the **2,000 most variable genes**. PCA gives every sample new coordinates (PC1, PC2…), ordered by the amount of variance they hold.
PCA is **unsupervised**: it knows nothing about donors, cell types or cytokines. To find out what a PC *means*, we **colour the samples by their metadata** and see which colouring lines up with the axes.
The **R²** panel quantifies it: the fraction of the variance of each PC that a variable explains (1 = completely).

**What to do:** the plot starts coloured by **donor**. Change the colouring (cell type, control vs treated, cytokine, number of cells, library size, a gene) and the PCs on the axes."""

EN["ex9_try"] = """#### ❓ Questions
1. With the default colouring (donor), do donors separate on PC1 and PC2? Colour by **cell type**: what do PC1 and PC2 separate, and how much variance do they hold?
2. Which PC best separates the donors? (try different PCs on the axes: the R² panel tells you how much each one has to do with the donors; then colour by donor.) And which PC has the largest R² with the cytokines?
3. Do cytokines have no effect on the cells, then? What does PC2's R² with the *library size* mean: should we worry?

<details>
<summary>💡 Answers</summary>

1. No: donors are mixed. Coloured by cell type, **PC1 (66 % of the variance)** separates the **monocytes** from the lymphocytes and **PC2 (17 %)** separates the **NK cells** from the T cells; CD4 and CD8 T cells are the closest pair. The R² of cell type with PC1 and PC2 is 0.99 and 0.98.
2. The donors separate on **PC6 (R² = 0.92)**, PC3 (0.72) and PC8 (0.65): try PC3 vs PC6. The cytokines have a modest R² on any single PC (at most 0.13, on PC3 and PC10): they matter, but on small PCs.
3. **No!** PCA shows the **largest** differences first, and the difference between a monocyte and a T cell is much bigger than the effect of any single cytokine. Cytokine effects live in later PCs or **inside each cell type** (next exercise). The R² of PC2 with the library size (0.43) is the cell type in disguise: NK cells simply have smaller libraries. After the normalization of Part C, depth is not driving the PCs.

</details>"""

EN["ex10"] = """## Exercise 10 · PCA within one cell type: cytokine responses
⏱ *10 min*

The cell-type differences hide almost everything else. A common remedy is to analyse **one cell type at a time**: we select its samples, choose its own most variable genes and run PCA again.
Within a cell type the donors and the cytokines are what varies.

**What to do:** choose a **cell type** (the default is monocytes, the most responsive) and the **PCs**. Colour by **cytokine** and tick the ones to highlight (PBS is always black), or colour by donor or by a gene.
The bar plots on the right show the **genes with the largest weights** (loadings) on each axis: they tell you what the axis *means*."""

EN["ex10_try"] = """#### ❓ Questions
1. In the default plot (PC1 vs PC4 in monocytes), which cytokines move away from the main cloud, and in which direction? Look at the genes behind each axis.
2. Switch the axes to PC1 vs PC2 and colour by **donor**. What dominates the first PCs within a cell type? Compare the R² of donor and cytokine in the table.
3. Colour by the expression of `UTY` (a Y-chromosome gene) and try several PCs: which PC separates male from female donors in monocytes? (the paper reports 6 males and 6 females)

<details>
<summary>💡 Answers</summary>

1. **IFN-β** moves to the right, along PC1, where the top genes include interferon-induced genes (`ACOD1`, `CCL8`, `CCL7`, `IDO1`); **IL-4 and IL-13** go up along PC4, whose top genes include `ALOX15` and `WNT5B` (IL-4 targets). **IL-4 and IL-13 land on the same place**: keep this in mind for Exercise 12.
2. **Donors**: the baseline differences between people are the largest variation within a cell type (R² of donor 0.72 on PC1 and 0.81 on PC2, vs 0.18 and 0.13 for cytokine). The cytokines show up more clearly in later PCs (R² 0.32 on PC3, 0.38 on PC4): the individual effects are smaller than the differences between people. That is why every cytokine will be compared with the PBS of the **same donor**.
3. **PC5** (correlation 0.87 between `UTY` and PC5). The sex of the donor is one of the largest sources of variation in monocytes, even if it has nothing to do with the cytokines. (In the other cell types it appears on another PC: PC2 in CD4 T cells, PC3 in CD8 T and NK cells.)

</details>"""

# ──────────────────────────────────────────────────────────────────────────────────────────────────
EN["partE"] = """---
# Part E · Differential expression: which genes change?

PCA showed that cytokines move the samples in particular directions. **Differential expression (DE)** asks the question **gene by gene**:

> *Is the expression of gene g with the cytokine different from its expression with PBS, by more than we would expect by chance?*

The results you will explore were **computed beforehand** with **DESeq2** (nothing is run now): for each of the 4 cell types and each of the 87 cytokines, the comparison *cytokine vs PBS* in the 12 donors.

> 🧬 **Only the immune genes.** To keep the analysis light and the plots readable, DE was run on **~500 immune-system genes** (506 expressed ones): the union of the GO sets *cytokine-mediated signaling pathway* and *T cell activation* and the Hallmark *interferon α / γ response* sets. Each comparison tests ~450 of them.
> **What do we lose?** (1) The list was built *from genes that respond to cytokines*, so the **fraction of DEGs is much higher** than in a whole-transcriptome analysis; do not compare percentages. (2) **Genes outside the list cannot be found**, even classic ones (`ACTB`, `GBP1`, `ALOX15`, `NKG7`, `CCL17`… are not in it): the tools tell you when a gene was not analysed.

<details>
<summary>🔬 How were the results computed? (optional)</summary>

* DESeq2 works on the **raw counts** and normalizes them itself (a robust version of CPM).
* The **design** `~ donor + cytokine` is **paired**: every donor is compared with itself, so differences between people cancel out.
* For each gene it fits a **negative binomial** model (counts are not normally distributed), and a **Wald test** asks whether the effect of the cytokine is larger than the noise.
* Because ~450 genes are tested at once, the p-values are **adjusted for multiple testing** (Benjamini–Hochberg, controlling the False Discovery Rate): the **adjusted p-value** is the one to use.
* Genes with too few counts (fewer than 10 counts in 12 of the 24 samples) are not tested.

</details>

Two numbers describe every gene:

| Number | Meaning |
|---|---|
| **log₂ fold change** | the **size** of the change: +1 = doubles (×2), +2 = ×4, −1 = halves |
| **adjusted p-value** | the **evidence**: how surprising the change would be if the cytokine did nothing (small = convincing) |

A gene is called **differentially expressed (DEG)** when it passes **both** thresholds.
"""

EN["ex11"] = """## Exercise 11 · Volcano plots: the genes that respond
⏱ *12 min*

A **volcano plot** shows every gene at once: *x* = log₂ fold change (effect size), *y* = −log₁₀ adjusted p-value (evidence; higher = more convincing).
**Red** = up-regulated DEGs, **blue** = down-regulated DEGs, **grey** = not significant. Genes at the top and far from the centre are the most robust and strongest.

**What to do:** the default is **IFN-γ in monocytes**, with the genes `CXCL9`, `CXCL10`, `IDO1`, `STAT1` and `IRF1` highlighted. Change the **cell type** and the **cytokine**, the **thresholds** and the genes to highlight
(type a symbol, or pick a ready-made list; only the ~500 immune genes can be found); the tables show their fold change, adjusted p-value and mean expression."""

EN["ex11_try"] = """#### ❓ Questions
1. How many genes go up and down with IFN-γ in monocytes? Look at the five highlighted genes: which are DEGs, and how many times higher than in PBS? Is any of them significant but *not* a DEG?
2. Set the fold-change threshold to 0, 1 and 2 (adjusted p-value 0.05), and then change the adjusted p-value threshold from 0.2 to 0.001 (fold change 0.5). Which threshold matters more? What does the fold-change threshold buy us?
3. Change the cytokine to **IL-4** and highlight `CCL26`, `CD1C`, `STAT1`, `IDO1`. What happens to the interferon genes?
4. Choose **NK cells** and **IL-2**. What do the top up-regulated genes `IL2RA`, `CISH` and `SOCS2` have in common?

<details>
<summary>💡 Answers</summary>

1. **96 up and 95 down** (of 452 genes tested). `CXCL9` (×10.5), `IDO1` (×5.5) and `CXCL10` (×3.7) are strongly up: the classic **IFN-γ response** (T-cell recruitment and tryptophan catabolism). `STAT1` is **not** a DEG (×1.1, adjusted p = 0.34), although it is the transcription factor of that response: it is already abundant in monocytes. `IRF1` is **significant** (adjusted p ≈ 10⁻⁵) but its change is small (×1.4, log₂FC 0.49 < 0.5): significant is not the same as large.
2. Fold-change threshold 0 / 1 / 2 → **311 / 70 / 14** DEGs; adjusted p-value 0.2 / 0.001 → 197 / 161 DEGs (with fold change 0.5; 191 with 0.05). The **fold-change threshold matters much more**: it removes statistically significant but biologically small changes, keeping those that are both reliable and large.
3. IL-4 strongly up-regulates its own targets (`CCL26` ×207, `CD1C` ×21) and **down-regulates** the interferon genes (`IDO1` ×0.1, `STAT1` ×0.4). IL-4 and IFN-γ push monocytes in opposite directions.
4. They are **IL-2 signalling genes**: `IL2RA` is the high-affinity IL-2 receptor α chain and `CISH` and `SOCS2` are STAT5 targets that feed back to switch the signal off. NK cells respond to IL-2 through STAT5.

</details>
"""

EN["ex12"] = """## Exercise 12 · Comparing two responses: similar or different?
⏱ *10 min*

Do two cytokines switch on the same genes? For every gene we have a log₂ fold change in **comparison A** and in **comparison B** (each = a cytokine vs PBS in one cell type). In a scatter plot of A against B,
two similar responses fall along the **diagonal**, and the **correlation coefficient** summarizes it: **r = 1** same response, **r = 0** unrelated, **r = −1** opposite.
Colours show where the gene is significant: **blue** only in A, **orange** only in B, **violet** in both, grey in neither. The right panel ranks all the other cytokines by how much they resemble A (in the same cell type).

**What to do:** the default compares **IL-4 and IL-13 in monocytes**. Change the cytokine and the cell type of A and B (they can be in different cell types!), load one of the **examples**, choose which genes are used for the correlation and **highlight genes**.

> ⚠️ Both comparisons use the **same PBS control**: noise in the PBS samples ends up in both and can create or inflate a correlation. Be careful with small values of r (≈ 0.2 or less)."""

EN["ex12_try"] = """#### ❓ Questions
1. How similar are the responses to IL-4 and IL-13 in monocytes? Now change both comparisons to **CD4 T cells**. What happens, and why? (remember `IL13RA1` from Exercise 7)
2. Load the examples **IFN-β vs IFN-γ in monocytes** and **IFN-γ vs IL-4 in monocytes**. Compare the correlations and use the right-hand panel: which cytokines respond most like IFN-γ?
3. Load **TNF-α in monocytes vs in CD4 T cells** and **IL-6 in monocytes vs in CD4 T cells**: which stimulus triggers a more similar response in the two cell types? And **TNF-α vs IL-6** in monocytes?
4. Change *genes used for the correlation* to "genes significant in A or B". When does r change a lot, and when not?

<details>
<summary>💡 Answers</summary>

1. In monocytes the two responses are **almost identical (r = 0.97)**, and IL-13 is the cytokine most similar to IL-4 out of 86. In CD4 T cells **r falls to 0.41** (IL-13 is only 3rd). IL-4 and IL-13 share a receptor (IL-4Rα/IL-13Rα1), but T cells hardly express `IL13RA1` (~2 CPM vs ~256 in monocytes), so they respond to IL-4 and not to IL-13. The **same pair of stimuli can look identical or very different depending on the cell type**.
2. **IFN-β vs IFN-γ: r = 0.78** (both interferons activate STAT1 and the interferon-stimulated genes); **IFN-γ vs IL-4: r = −0.35**: opposite programmes (what one raises, the other lowers). The ranking for IFN-γ in monocytes starts with IL-12 (0.87), IFN-β and IFN-ω (0.78) and IL-15 (0.74), and ends with IL-17E and TNF-α (−0.46, −0.45).
3. **TNF-α** gives a similar response in the two cell types (**r = 0.70**), **IL-6** a more cell-type-specific one (**r = 0.47**). **TNF-α vs IL-6 in monocytes: r = 0.22**, a weak relation. Always look at the data before assuming that two cytokines share a program.
4. In monocytes, where both cytokines change many genes, r hardly moves (IL-4 vs IL-13: 0.97 → 0.97). Where responses are weaker, restricting to genes that respond (CD4 T cells: 0.41 → 0.51; TNF-α vs IL-6 in monocytes: 0.22 → 0.28) removes noise and r **increases**. But it also *selects* genes that change, so the comparison is no longer free of bias: report which genes were used. Remember too that the whole analysis already uses a list of immune genes, enriched in cytokine-responsive genes: correlations are higher than they would be genome-wide.

</details>
"""

# ──────────────────────────────────────────────────────────────────────────────────────────────────
EN["wrapup"] = """---
# 🏁 Wrap-up

You went from a single variable to the comparison of transcriptional programs:

| Exercise | The idea to keep |
|:--:|---|
| 1–2 | A variable is a name for a value; a plot type answers a question. |
| 3 | AnnData keeps the numbers (`X`), the samples (`obs`) and the genes (`var`) **aligned**; `adata[samples, genes]` selects a subset. |
| 4 | Raw counts mix biology with the **amount of material** (library size). |
| 5 | **CPM** divides by the library size so that samples become comparable. |
| 6 | The **logarithm** compresses the scale, turns ratios into differences and stops a few genes from dominating. |
| 7 | **Heatmaps** show patterns; they only make sense *after* normalizing. |
| 8–10 | **PCA** is a rotation that finds the directions of largest variance. Colouring by metadata tells what each PC means: cell type ≫ donor ≫ cytokine. Look **inside a cell type** to see cytokine effects. |
| 11 | A **DEG** needs evidence (**adjusted p-value**) *and* a relevant effect (**fold change**). The paired design compares each donor with itself. |
| 12 | The **correlation of fold changes** says whether two stimuli trigger similar programs, and it **depends on the cell type**. |

### Where to go from here
* **Functional enrichment** of the lists of DEGs (which pathways?): `gseapy`, `decoupler`.
* **Single-cell analysis** of the original data: `scanpy`; best practices in Heumos *et al.* (2023).
* The **Human Cytokine Dictionary** itself and its `huCIRA` package (Oesinghaus *et al.*, 2025) to interpret your own data with cytokine-response signatures.

### 📝 Quick quiz
Your instructor will show the **QR code** of the quiz (3 questions, 1 minute).

### References
* Oesinghaus L. *et al.* (2025) A single-cell cytokine dictionary of human peripheral blood. *bioRxiv*, doi:10.64898/2025.12.12.693897
* Virshup I. *et al.* (2024) anndata: Access and store annotated data matrices. *J. Open Source Softw.* 9:4371
* Heumos L. *et al.* (2023) Best practices for single-cell analysis across modalities. *Nat. Rev. Genet.* 24:550–572
* Love M.I., Huber W., Anders S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.* 15:550
* Muzellec B. *et al.* (2023) PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* 39:btad547
* Benjamini Y., Hochberg Y. (1995) Controlling the false discovery rate. *J. R. Stat. Soc. B* 57:289–300
"""

TEXT = {"en": EN}
