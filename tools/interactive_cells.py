"""Code of the cells of the interactive notebook (language independent).

The helper library itself lives in tools/lib/*.py and is pasted into the hidden setup cell; the exercise cells only call one function."""
from libsource import library_source
from translate_lib import translate_source
from translations_es import ES_DISPLAY, ES_STRINGS

INSTALL = '''# ── Install whatever is missing (Google Colab already has most packages) ──────────────────────────
import importlib.util, shutil, subprocess, sys

REQUIRED = {                       # import name -> PyPI package name
    "numpy": "numpy", "pandas": "pandas", "scipy": "scipy", "matplotlib": "matplotlib", "anndata": "anndata",
    "ipywidgets": "ipywidgets", "jinja2": "jinja2", "gdown": "gdown", "adjustText": "adjustText",
}
missing = [pkg for mod, pkg in REQUIRED.items() if importlib.util.find_spec(mod) is None]
if missing:
    print(f"📦 Installing: {', '.join(missing)} …")
    if shutil.which("uv"):
        cmd = ["uv", "pip", "install", "--python", sys.executable, *missing]
    else:
        cmd = [sys.executable, "-m", "pip", "install", "-q", *missing]
    subprocess.check_call(cmd)
    importlib.invalidate_caches()
'''

SETUP_TITLE = '#@title ▶ Setup: run this cell first (about 2 minutes) { display-mode: "form" }'


SETUP_TITLE_ES = '#@title ▶ Preparación: ejecuta esta celda primero (unos 2 minutos) { display-mode: "form" }'
MISSING = {"en": "⚠️ Run the Setup cell at the top of the notebook first (▶), wait for ✅ Ready, and then run this cell again.",
           "es": "⚠️ Ejecuta primero la celda de Preparación al principio del cuaderno (▶), espera a ✅ Listo y vuelve a ejecutar esta celda."}


def lib_source(lang="en") -> str:
    src = library_source()
    if lang == "es":
        src, used = translate_source(src, ES_STRINGS)
        unused = [k for k in ES_STRINGS if k not in used]
        if unused:
            print(f"ℹ️  {len(unused)} Spanish strings not found in the library (harmless, but check typos): {[u[:30] for u in unused[:8]]}")
        src += f"\nDISPLAY.update({ES_DISPLAY!r})\n"
    return src


def setup_cell_source(lang="en") -> str:
    install = INSTALL if lang == "en" else INSTALL.replace("Install whatever is missing (Google Colab already has most packages)", "Instala lo que falte (Google Colab ya trae casi todo)").replace("📦 Installing:", "📦 Instalando:")
    return "\n".join([SETUP_TITLE if lang == "en" else SETUP_TITLE_ES, install, "", lib_source(lang), "", "# ── Load the data and pre-compute everything ──────────────────────────────────────────────────────",
                      "init_practical()", ""])


def exercise_cell_source(title: str, function: str, lang="en") -> str:
    """One-line exercise cell, with a friendly message if the setup cell was not run."""
    return (f'#@title ▶ {title} {{ display-mode: "form" }}\n'
            f"try:\n    {function}()\n"
            f"except NameError:\n    print({MISSING[lang]!r})\n")


# (markdown key, cell title, function, "try this" key)
EXERCISES = [
    ("ex1", "Exercise 1 · Variables", "exercise_variables", "ex1_try"),
    ("ex2", "Exercise 2 · Basic plots", "exercise_plotting", "ex2_try"),
    ("ex3", "Exercise 3 · The AnnData object", "exercise_anndata", "ex3_try"),
    ("ex4", "Exercise 4 · Boxplots of raw counts", "exercise_boxplots", "ex4_try"),
    ("ex5", "Exercise 5 · Raw counts vs CPM", "exercise_normalization", "ex5_try"),
    ("ex6", "Exercise 6 · The logarithm", "exercise_logtransform", "ex6_try"),
    ("ex7", "Exercise 7 · Heatmaps", "exercise_heatmaps", "ex7_try"),
    ("ex8", "Exercise 8 · PCA, the idea", "exercise_pca_geometry", "ex8_try"),
    ("ex9", "Exercise 9 · PCA of all samples", "exercise_pca_all", "ex9_try"),
    ("ex10", "Exercise 10 · PCA within a cell type", "exercise_pca_celltype", "ex10_try"),
    ("ex11", "Exercise 11 · Volcano plots", "exercise_volcano", "ex11_try"),
    ("ex12", "Exercise 12 · Comparing two responses", "exercise_compare", "ex12_try"),
]
TITLES_ES = ["Ejercicio 1 · Variables", "Ejercicio 2 · Gráficos básicos", "Ejercicio 3 · El objeto AnnData", "Ejercicio 4 · Diagramas de cajas de conteos crudos",
             "Ejercicio 5 · Conteos crudos vs CPM", "Ejercicio 6 · El logaritmo", "Ejercicio 7 · Mapas de calor", "Ejercicio 8 · PCA, la idea",
             "Ejercicio 9 · PCA de todas las muestras", "Ejercicio 10 · PCA dentro de un tipo celular", "Ejercicio 11 · Volcano plots", "Ejercicio 12 · Comparar dos respuestas"]
PARTS = {"ex1": "partA", "ex3": "partB", "ex5": "partC", "ex7": "partD", "ex11": "partE"}      # a part header goes before these exercises
