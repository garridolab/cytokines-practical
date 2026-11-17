"""Narrativa en español del cuaderno interactivo (mismas claves que interactive_text.EN). Cifras comprobadas con los datos."""
ES = {}

ES["title"] = """# 🧬 Diccionario humano de citocinas: una práctica interactiva de transcriptómica

{badge}

**Basado en:** Oesinghaus *et al.* (2025), *A single-cell cytokine dictionary of human peripheral blood*, bioRxiv
([doi:10.64898/2025.12.12.693897](https://doi.org/10.64898/2025.12.12.693897)).

En el artículo, ~9,7 millones de células sanguíneas humanas (PBMC) de **12 donantes** se estimularon con **90 citocinas** (o con un control de PBS)
y se midieron por RNA-seq de célula única. Aquí trabajas con una versión **pseudobulk** de los datos: para cada *donante × citocina × tipo celular*
se **sumaron** los conteos de todas las células en un único perfil de expresión tipo bulk. Usamos **cuatro tipos celulares** (células T CD4, células T CD8, células NK, monocitos CD14),
**87 citocinas** y el **control PBS**.

### 🎛 No vas a escribir código
Esta práctica es **interactiva**. Cada ejercicio te da unos pocos **controles** (menús, deslizadores, cajas) y una **figura que se redibuja sola** cada vez que cambias uno.
Tu trabajo es **jugar, observar e interpretar**: ¿qué cambia? ¿por qué?

| Parte | Ejercicios | Tiempo |
|:--:|---|:--:|
| **A · Calentamiento en Python** | **1** Variables · **2** Gráficos básicos | 12 min |
| **B · Conoce los datos** | **3** El objeto AnnData · **4** Distribuciones con diagramas de cajas | 18 min |
| **C · Normalización** | **5** Conteos crudos vs CPM · **6** El logaritmo | 18 min |
| **D · Patrones en los datos** | **7** Mapas de calor · **8** PCA: la idea · **9** PCA: todas las muestras · **10** PCA: un tipo celular | 38 min |
| **E · Expresión diferencial** | **11** Volcano plots · **12** Comparar dos respuestas | 22 min |
| | Cierre y cuestionario | 7 min |

### Cómo usar este cuaderno

| | |
|:--:|---|
| ▶ | **Ejecuta una celda** con el botón ▶ de su izquierda (o haz clic y pulsa **Shift + Enter**). Ejecuta primero la **celda de preparación** y después cada celda de ejercicio *una vez*: aparecen los controles. |
| 🎛 | **Cambia un control y espera un segundo**: figuras y tablas se actualizan solas. No hace falta volver a ejecutar la celda. |
| 🖱 | En una lista, **Ctrl/⌘ + clic** selecciona varios elementos. En las cajas de texto, pulsa **Enter** al terminar de escribir. |
| ❓ | Cada ejercicio termina con **preguntas**. Piensa (o coméntalo con tu compañero) *antes* de abrir la 💡 respuesta. |
| 🔁 | Si algo se ve raro, vuelve a ejecutar la celda del ejercicio. Si el cuaderno se reinició, ejecuta de nuevo la celda de preparación. |

> 💡 En Google Colab el código de las celdas grises está oculto a propósito. Si tienes curiosidad, haz doble clic en el título de la celda para verlo.
> Los números usan punto decimal (0.93) y coma para los millares (1,353), igual que las figuras.
"""

ES["setup_md"] = """---
# ⚙️ Preparación: ejecuta esta celda primero

Instala lo que falte, descarga los datos (≈ 500 MB, alrededor de un minuto en Colab), los normaliza, calcula los PCA y carga los resultados de expresión diferencial
calculados de antemano con DESeq2. **Solo la ejecutas una vez** (y otra vez si el cuaderno se reinicia). Espera a que aparezca **✅ Listo**."""

ES["partA"] = """---
# Parte A · Calentamiento en Python
Dos ejercicios cortos que no necesitan datos: qué es una **variable** y cómo **mirar números** con gráficos."""

ES["ex1"] = """## Ejercicio 1 · Variables: ponerle nombre a los valores
⏱ *5 min*

Una **variable** es un **nombre** que guarda un **valor**, como un tubo etiquetado en el congelador. Se crea con `=`: el nombre a la izquierda, el valor a la derecha.
Python necesita saber *de qué tipo* es el valor: **texto** (`str`, siempre entre comillas), **números enteros** (`int`), **números decimales** (`float`), **Verdadero/Falso** (`bool`) o una **lista** de valores.
`print(...)` muestra un valor en pantalla, y las variables se pueden combinar para crear otras nuevas.

**Qué hacer:** cambia los controles y observa cómo cambian el código Python, lo que se imprime y la tabla de variables."""

ES["ex1_try"] = """#### ❓ Preguntas
1. Cambia `message` por tu nombre. ¿Qué líneas del panel de código y de la salida cambian?
2. Mueve el deslizador `n_donors`. ¿Qué líneas impresas cambian y cuáles no? ¿Por qué?
3. En la tabla de variables, ¿de qué tipo es `library_size`? ¿Por qué no es un número entero?

<details>
<summary>💡 Respuestas</summary>

1. La línea `message = '...'` y la primera línea impresa. Los textos se escriben **entre comillas**; la variable guarda lo que escribas.
2. Solo las líneas que *usan* `n_donors`: `Donantes: …` y `… se añadió a … donantes`. Cambiar una variable solo afecta a los sitios donde aparece su nombre.
3. Es un `float` (decimal): `library_size = n_cells * umis_per_cell` multiplica un entero por un decimal y el resultado es decimal. Si ambos fueran enteros (`int`), el resultado sería `int`. (`library_size`, el número total de conteos de una muestra, volverá a aparecer en el Ejercicio 3.)

</details>"""

ES["ex2"] = """## Ejercicio 2 · Mirar números: gráficos básicos
⏱ *7 min*

Un gráfico convierte muchos números en una imagen. Cada tipo de gráfico responde a una pregunta distinta:

| Gráfico | Responde a la pregunta… |
|---|---|
| 📊 **Histograma** | *¿Cómo se distribuyen los valores?* (cuántos donantes tienen cada valor) |
| 📦 **Diagrama de cajas** | *¿Dónde está el centro y cuánto se dispersan los valores?* Línea = mediana, caja = el 50 % central, bigotes = el rango |
| 📶 **Barras** | *¿Cuál es la media de cada grupo?* (¡oculta los valores individuales!) |
| 🔵 **Dispersión** | *¿Están relacionadas dos medidas?* Un punto por donante |
| 📈 **Líneas** | *¿Cómo cambia un valor con el tiempo?* |

**Qué hacer:** los datos son **simulados**: un gen medido en varios donantes, sin (control) y con (tratado) una citocina. Elige un tipo de gráfico y juega con los números."""

ES["ex2_try"] = """#### ❓ Preguntas
1. En el **histograma**, pon el efecto de la citocina a 0 y luego a +2. ¿Qué se mueve? ¿Qué pasa si además aumentas la variabilidad entre donantes (DE)?
2. Usa pocos donantes (5–10) y pulsa **🎲 nuevo sorteo aleatorio** varias veces (histograma o cajas). ¿Los dos grupos se ven siempre igual? ¿Qué te dice esto de los experimentos con pocos donantes?
3. En el **diagrama de dispersión**, usa un efecto de +1 y compara un parecido de 0.1 con otro de 0.95. ¿Cuándo se ve mejor el efecto? ¿Por qué ayudaría comparar cada donante *consigo mismo*?

<details>
<summary>💡 Respuestas</summary>

1. Con efecto 0 los dos histogramas se solapan por completo; con +2 el tratado se desplaza 2 unidades a la derecha (2 unidades de log₂ = **4 veces** más expresión). Una DE mayor ensancha ambos histogramas, que se solapan más aunque el efecto sea el mismo: el efecto es más difícil de ver cuando los donantes difieren mucho.
2. No: con pocos donantes el azar cambia las medias y las formas, y puedes ver "diferencias" (o perder reales) por mala suerte. Con más donantes las imágenes son estables. Por eso los experimentos necesitan réplicas.
3. Con parecido 0.95 todos los puntos quedan en una banda estrecha **por encima** de la línea discontinua (`y = x`): *cada* donante sube más o menos lo mismo, aunque partan de niveles muy distintos. Con 0.1 la nube es redonda y el desplazamiento queda tapado por las diferencias entre donantes. Comparar cada donante consigo mismo (un diseño **pareado**, Ejercicios 11–12) elimina las diferencias entre personas.

</details>"""

ES["partB"] = """---
# Parte B · Conoce los datos"""

ES["ex3"] = """## Ejercicio 3 · El objeto AnnData: muestras, genes y sus metadatos
⏱ *10 min*

Nuestros datos viven en un objeto **AnnData** llamado `adata`, que mantiene varias tablas **alineadas**:

| Componente | Qué contiene | Forma |
|---|---|---|
| `X` | los **números**: conteos crudos (cuántas moléculas de mRNA de cada gen se detectaron) | muestras × genes |
| `obs` | información de las **muestras** (filas de `X`): donante, citocina, tipo celular… | una fila por muestra |
| `var` | información de los **genes** (columnas de `X`) | una fila por gen |
| `layers` | otras versiones de `X` con la misma forma (aquí, versiones *normalizadas* que construiremos en los Ejercicios 5–6) | muestras × genes |

El objeto completo tiene **4,224 muestras × 18,120 genes**: 12 donantes × 4 tipos celulares × 88 condiciones (87 citocinas + PBS).
Una **muestra** es un pseudobulk: todas las células de un donante, un tipo celular y una condición, sumadas. Su nombre es `donante|tipo celular|citocina`.
Seleccionar funciona como en una tabla: `adata[muestras, genes]`, **primero muestras, después genes**.

**Qué hacer:** elige muestras y genes (por defecto: los 4 controles PBS del Donor 1, más sus monocitos tratados con IFN-γ). El **subconjunto** que construyes aparece en las tres pestañas, y su lugar en la matriz completa, a la derecha."""

ES["ex3_try"] = """#### ❓ Preguntas
1. En el subconjunto por defecto, ¿qué muestra tiene el mayor conteo de `LYZ`, de `CD3E` y de `NKG7`? ¿Concuerda con lo que sabes de estos genes?
2. Abre la pestaña **obs**: ¿qué te dicen `psbulk_n_cells` y `library_size`? Compara las dos muestras de monocitos (PBS e IFN-γ).
3. `STAT1` se activa con IFN-γ. Compara su conteo en las dos muestras de monocitos. ¿Es mayor con IFN-γ? ¿Es lo que esperabas?

<details>
<summary>💡 Respuestas</summary>

1. `LYZ` es máximo en los **monocitos** (911,945 conteos), `CD3E` en las **células T CD4** (20,361; las T CD8 son segundas con 3,100) y `NKG7` en las **células NK** (2,392; luego las T CD8, también citotóxicas): son marcadores de tipo celular.
2. `psbulk_n_cells` es el **número de células que se sumaron** en la muestra y `library_size` el **número total de conteos** (de todos los genes). La muestra PBS de monocitos contiene 17,373 células y 147.9 millones de conteos; la de IFN-γ solo 2,556 células y 21.4 millones: unas **7 veces menos material**. (El control PBS es la suma de 6 pocillos.)
3. **No: es ~9 veces *menor* con IFN-γ** (73,853 frente a 639,289), aunque IFN-γ induce `STAT1`. El conteo crudo mezcla **biología** (¿está el gen activo?) con la **cantidad de material** (cuántas células, qué profundidad). Es el problema que resolvemos en la Parte C.

</details>"""

ES["ex4"] = """## Ejercicio 4 · Distribuciones de conteos crudos: diagramas de cajas
⏱ *8 min*

Un **diagrama de cajas** resume muchos valores: la **línea** es la mediana, la **caja** contiene el 50 % central de los valores y los **bigotes** muestran el rango. Cada **punto** es una **muestra pseudobulk**
(un donante con una citocina en un tipo celular); la *caja* es, por tanto, la distribución del gen en todas las muestras de un grupo.
El panel derecho muestra el **tamaño de librería** (conteos totales) de las mismas muestras.

**Qué hacer:** por defecto verás cuatro genes con distintos papeles en los cuatro tipos celulares (las 88 condiciones, 1,056 muestras por tipo celular). Cambia los genes, los grupos o los tipos celulares.
Prueba "control (PBS) vs tratado" para comparar los dos tipos de muestras."""

ES["ex4_try"] = """#### ❓ Preguntas
1. En el gráfico por defecto, ¿cuántos órdenes de magnitud separan el gen más y el menos abundante en monocitos? ¿Se comportan `LYZ`, `CD3E` y `NKG7` como marcadores de tipo celular?
2. `ACTB` es un gen constitutivo (siempre activo), pero sus conteos difieren entre tipos celulares. Compara con el panel derecho: ¿hay relación?
3. Agrupa las cajas por **control (PBS) vs tratado**. ¿Qué ves para *todos* los genes? ¿Por qué?

<details>
<summary>💡 Respuestas</summary>

1. En monocitos la mediana va de 4 (`NKG7`) a 25,725 (`STAT1`): **casi 4 órdenes de magnitud**. `LYZ` es máximo en monocitos (mediana 30,058 frente a ≤ 297 en el resto), `CD3E` en células T CD4 y CD8 (1,262 y 663) y `NKG7` en NK (134; luego T CD8, 65).
2. Sí: el **tamaño de librería** difiere entre tipos celulares (medianas: NK 1.2 M, T CD8 3.8 M, monocitos 7.6 M, T CD4 8.7 M) y los conteos de `ACTB` lo siguen en buena parte (NK el menor, luego T CD8). Los monocitos son células más grandes con más mRNA, así que tienen aún más `ACTB` por conteo. Los conteos dependen de cuánto material aporta cada tipo celular.
3. Los **controles PBS son más altos en todos los genes**, incluido `ACTB`: su tamaño de librería es ~5 veces mayor (mediana 23 M frente a 4.4 M) porque cada muestra PBS es la suma de 6 pocillos. Esa diferencia no tiene que ver con la biología de cada gen, sino con la cantidad de material, que corregimos en el siguiente ejercicio.

</details>"""

ES["partC"] = """---
# Parte C · Normalización: hacer comparables las muestras"""

ES["ex5"] = """## Ejercicio 5 · Conteos crudos vs CPM: ¿por qué normalizamos?
⏱ *10 min*

Un conteo crudo depende de **dos cosas**: cuánto se expresa el gen (biología) y cuánto material tenemos (cuántas células, qué profundidad de secuenciación; técnica).
El **tamaño de librería** es el número total de conteos de una muestra. Para quitar el efecto de la profundidad expresamos cada gen como fracción de su muestra:

**CPM** (conteos por millón) = conteo del gen ÷ tamaño de librería de la muestra × 1,000,000

**Qué hacer:**
* **Ejemplo de juguete** (simulado): el mismo donante con una muestra profunda (PBS) y otra superficial (con una citocina). Solo el gen B cambia de verdad. Mueve los deslizadores.
* **Datos reales, PBS vs citocina**: una línea por donante. Por defecto, el gen constitutivo `ACTB` en monocitos con IFN-γ. Cambia el gen, el tipo celular o la citocina.
* **Datos reales, todas las muestras**: ¿sigue el conteo del gen al tamaño de librería?"""

ES["ex5_try"] = """#### ❓ Preguntas
1. **Ejemplo de juguete**: con una profundidad del 20 % y un cambio real de ×4, ¿qué dicen los conteos crudos del gen B? ¿Y el CPM? ¿Y con profundidad del 100 %?
2. **Datos reales**, `ACTB` en monocitos (PBS vs IFN-γ): ¿cómo cambian el tamaño de librería, el conteo crudo de `ACTB` y su CPM? ¿Qué esperas de un gen constitutivo?
3. Cambia el gen a `IDO1` y luego a `STAT1` (ambos inducidos por IFN-γ). ¿Qué dicen conteos crudos y CPM? Después mira la vista **todas las muestras**: ¿cómo cambia la correlación con el tamaño de librería de crudos a CPM?

<details>
<summary>💡 Respuestas</summary>

1. **Conteos crudos: el gen B parece ligeramente *a la baja* (×0.78)** y todos los demás genes parecen ~5 veces más bajos, solo porque la muestra 2 tiene 5 veces menos profundidad. **CPM: los genes A, C y D no cambian (×0.97) y el gen B es ×3.9**, cerca del ×4 real. Con profundidad del 100 % ambas representaciones coinciden: normalizar solo importa cuando las profundidades difieren.
2. El tamaño de librería es **6 veces menor** con IFN-γ (×0.16) y el conteo crudo de `ACTB` lo sigue (×0.16); tras el CPM `ACTB` **no cambia (×0.95)**, como se espera de un gen constitutivo. En la vista de todas las muestras la correlación de `ACTB` con el tamaño de librería es **r = +0.84 en crudos y −0.25 en CPM**: el efecto de la profundidad desaparece.
3. `IDO1`: los conteos crudos dicen **sin cambio (×1.1)**, el CPM dice **×6.7 al alza**. `STAT1`: los crudos dicen **×4.5 *a la baja***, el CPM **×1.35 al alza**. Sin normalizar perderíamos la inducción de `IDO1` y concluiríamos lo contrario para `STAT1`.

</details>"""

ES["ex6"] = """## Ejercicio 6 · El logaritmo: domar la escala
⏱ *8 min*

El CPM arregla las diferencias entre **muestras**, pero no entre **genes**: unos tienen <1 CPM y otros >10,000, los valores son muy **asimétricos** y la **varianza crece con la media**.
Muchos análisis (PCA, tests…) se guían por la varianza, así que unos pocos genes abundantes lo decidirían todo. El **logaritmo** comprime los valores grandes y convierte cocientes en diferencias
(con log₂: duplicar = +1, reducir a la mitad = −1). Como log(0) no existe, primero sumamos un pequeño **pseudoconteo** (normalmente 1): **log₂(CPM + 1)**.

**Qué hacer:** los paneles superiores muestran un gen en un tipo celular; los inferiores, **todos los genes** (varianza frente a media). Cambia la **base** del logaritmo y el **pseudoconteo**."""

ES["ex6_try"] = """#### ❓ Preguntas
1. Compara ① y ② para `CXCL10` en monocitos. ¿Qué le pasa a la forma de la distribución y a su varianza?
2. Compara ③ y ④: ¿qué parte de la varianza total acumulan los 20 genes más variables, antes y después del logaritmo? ¿Qué significaría eso para un PCA?
3. Cambia la base a 10 y a *e*. ¿Qué cambia y qué no? ¿Y con el pseudoconteo a 0 y a 10?

<details>
<summary>💡 Respuestas</summary>

1. En CPM la distribución es **asimétrica** (cola larga: media 1,297, mediana 634, asimetría 1.6) y su varianza enorme (2,350,366). Tras log₂ es mucho más simétrica (asimetría −0.3) y la varianza es 4.7.
2. En CPM **20 genes acumulan el 41 % de toda la varianza**; tras log₂(CPM + 1) solo el **3 %**. Sin logaritmo, un PCA lo decidirían unos pocos genes abundantes (y casi solo "vería" esos); con él contribuyen miles de genes.
3. La **forma no cambia**, solo los números del eje (log₂ x = 3.32 × log₁₀ x): las varianzas se multiplican por 0.48 (ln) o 0.09 (log₁₀). La base es una convención; usamos 2 porque una diferencia de 1 significa "el doble". Con pseudoconteo **0**, los ceros no se pueden transformar y se descartan. Con **10**, los valores bajos se aplastan. 1 es el compromiso habitual.

</details>"""

ES["partD"] = """---
# Parte D · Buscar patrones en los datos
A partir de ahora usamos los valores **normalizados**: log₂(CPM + 1)."""

ES["ex7"] = """## Ejercicio 7 · Mapas de calor: genes × muestras
⏱ *10 min*

Un **mapa de calor** muestra muchos genes (filas) en muchas muestras (columnas) a la vez. Para comparar patrones, cada gen (fila) se **escala entre las muestras** (z-score): **rojo** = por encima de la media de *ese* gen, **azul** = por debajo.
Las dos barras superiores muestran el **tipo celular** y el **tamaño de librería** de cada columna (oscuro = profundo).

**Qué hacer:** por defecto verás marcadores de tipo celular y genes de respuesta a citocinas en los cuatro tipos celulares con PBS y cuatro citocinas (12 donantes cada una), con los **conteos crudos** a la izquierda y los **valores normalizados** a la derecha.
Cambia de vista, elige otros **grupos de genes**, añade tus propios genes (p. ej. `TNF`, `IL1B`, `IL13RA1`) y selecciona otras **citocinas** (Ctrl/⌘ + clic)."""

ES["ex7_try"] = """#### ❓ Preguntas
1. Compara ANTES y DESPUÉS. ¿Qué columnas son las más brillantes en el mapa crudo, para casi todos los genes? ¿Puedes leer los tipos celulares en cada panel?
2. En el mapa normalizado, ¿qué citocina activa qué bloque de genes, y en qué tipo celular?
3. Selecciona **IFN-beta** y el grupo *Respuesta a interferón tipo I*. ¿En qué tipos celulares responden estos genes? Después añade `IL13RA1` e `IL4R` como genes extra: ¿qué tipo celular expresa `IL13RA1`? (tenlo presente para el Ejercicio 12)

<details>
<summary>💡 Respuestas</summary>

1. **Antes**, las columnas más brillantes son las **PBS** (oscuras en la barra de tamaño de librería), sea cual sea el gen: la imagen muestra *cuántos conteos tiene cada muestra* y los tipos celulares cuesta leerlos. **Después**, el brillo de PBS desaparece y cada tipo celular muestra sus marcadores: `LYZ`, `CD14`, `S100A9`, `FCN1` en monocitos; `CD3E`, `CD3D` en células T; `CD8A`, `CD8B` en T CD8; `NKG7`, `GNLY`, `KLRF1` en NK.
2. **IL-4** → `ALOX15`, `CCL17`, `FCER2` casi solo en **monocitos**. **IFN-γ** → `GBP1`, `IDO1`, `CXCL9`, `CXCL10`, máximo en **monocitos**. **IL-1β** → `IL6`, `CXCL8`, `CCL2` en **monocitos**. **IL-2** → `IL2RA`, `GZMB`, `CISH` en **NK y T**, los linfocitos que tienen el receptor de IL-2.
3. Los genes estimulados por interferón (`ISG15`, `IFIT1`, `MX1`…) suben **en los cuatro tipos celulares** con IFN-β (log₂ fold change entre 2 y 3), porque el receptor de interferón tipo I está en todos. `IL13RA1` (una subunidad del receptor de IL-13) se expresa en **monocitos** (~256 CPM) y prácticamente no en T y NK (~2 CPM).

</details>

> ⚠️ **Los z-score son para visualizar.** Ocultan el nivel absoluto de cada gen y exageran genes que apenas cambian. Nunca calcules fold changes a partir de ellos."""

ES["ex8"] = """## Ejercicio 8 · PCA: la idea, con solo dos genes
⏱ *10 min*

El **PCA** (análisis de componentes principales) resume muchos genes en unos pocos ejes nuevos. Para ver cómo, toma solo **dos genes**: cada muestra es un punto en un plano (*x* = gen 1, *y* = gen 2).
1. Dibuja un **eje** (una recta por el centro de la nube) y **proyecta** cada muestra sobre él (baja una perpendicular desde el punto hasta la recta).
2. La **dispersión (varianza)** de los puntos proyectados indica cuánta información conserva ese eje.
3. **Gira** el eje hasta que la varianza sea **máxima**: ese eje es la **PC1**. El perpendicular es la **PC2**.

El PCA es *solo un giro* de los ejes (no se pierde información); la **reducción de dimensionalidad** aparece cuando conservamos solo las primeras PC. Con miles de genes la receta es exactamente la misma.

**Qué hacer:** hay dos ejemplos fijos en monocitos: un par de genes **correlacionados** (`ISG15` y `MX1`, ambos activados por interferón) y un par **no correlacionados** (`IDO1` e `IL1B`, dos rutas distintas).
Gira el eje con el deslizador y observa los paneles ①②③. Marca *mostrar PC1 y PC2* o pulsa *ir a PC1*."""

ES["ex8_try"] = """#### ❓ Preguntas
1. **Par correlacionado**: ¿con qué ángulo conserva el eje más varianza? ¿Cuánta? ¿Cómo es el eje respecto a la nube de puntos?
2. **Par no correlacionado**: ¿cuánto cambia la varianza conservada con el ángulo? ¿Podría un eje resumir ambos genes?
3. En el panel ③, con el ángulo de la PC1, ¿cuál es la correlación entre las dos nuevas coordenadas? Si descartamos la PC2, ¿cuánta varianza perdemos en cada par?

<details>
<summary>💡 Respuestas</summary>

1. Con unos **36°** el eje conserva el **98 %** de la varianza total: corre a lo largo de la dirección larga de la nube alargada (los dos genes suben y bajan juntos, r = 0.95). Un número por muestra (la posición a lo largo de la PC1) puede sustituir a los dos genes.
2. Muy poco: entre el **45 % y el 55 %** sea cual sea el ángulo (el mejor eje, la PC1, conserva el 55 %). La nube no tiene dirección preferida, así que ningún eje resume ambos genes: perderíamos casi la mitad de la varianza.
3. **Alrededor de 0**: PC1 y PC2 están incorrelacionadas por construcción. Descartar la PC2 pierde solo el **2 %** de la varianza en el par correlacionado, pero el **45 %** en el no correlacionado. El PCA comprime bien *cuando los genes están correlacionados*: con miles de genes que responden juntos (todos los de interferón, todos los de monocito…) unas pocas PC resumen la mayor parte de la variación.

</details>"""

ES["ex9"] = """## Ejercicio 9 · PCA de todas las muestras: ¿qué determina la variación?
⏱ *8 min*

Cada una de las 4,224 muestras se describe con los **2,000 genes más variables**. El PCA da a cada muestra nuevas coordenadas (PC1, PC2…), ordenadas por la varianza que contienen.
El PCA es **no supervisado**: no sabe nada de donantes, tipos celulares o citocinas. Para saber qué *significa* una PC, **coloreamos las muestras por sus metadatos** y vemos qué coloración se alinea con los ejes.
El panel de **R²** lo cuantifica: la fracción de la varianza de cada PC que explica una variable (1 = por completo).

**Qué hacer:** el gráfico empieza coloreado por **donante**. Cambia la coloración (tipo celular, control vs tratado, citocina, número de células, tamaño de librería, un gen) y las PC de los ejes."""

ES["ex9_try"] = """#### ❓ Preguntas
1. Con la coloración por defecto (donante), ¿se separan los donantes en PC1 y PC2? Colorea por **tipo celular**: ¿qué separan PC1 y PC2 y cuánta varianza contienen?
2. ¿Qué PC separa mejor a los donantes? (prueba distintas PC en los ejes: el panel de R² indica cuánto tiene cada una que ver con los donantes; luego colorea por donante.) ¿Y qué PC tiene mayor R² con las citocinas?
3. ¿Entonces las citocinas no tienen efecto sobre las células? ¿Qué significa el R² de la PC2 con el *tamaño de librería*: debe preocuparnos?

<details>
<summary>💡 Respuestas</summary>

1. No: los donantes están mezclados. Coloreado por tipo celular, la **PC1 (66 % de la varianza)** separa los **monocitos** de los linfocitos y la **PC2 (17 %)** separa las **células NK** de las T; T CD4 y T CD8 son el par más cercano. El R² del tipo celular con PC1 y PC2 es 0.99 y 0.98.
2. Los donantes se separan en la **PC6 (R² = 0.92)**, la PC3 (0.72) y la PC8 (0.65): prueba PC3 vs PC6. Las citocinas tienen un R² modesto en cualquier PC (como mucho 0.13, en PC3 y PC10): importan, pero en PC pequeñas.
3. **¡No!** El PCA muestra primero las diferencias **mayores**, y la diferencia entre un monocito y una célula T es mucho mayor que el efecto de cualquier citocina. Los efectos de las citocinas están en PC posteriores o **dentro de cada tipo celular** (siguiente ejercicio). El R² de la PC2 con el tamaño de librería (0.43) es el tipo celular disfrazado: las NK simplemente tienen librerías más pequeñas. Tras la normalización de la Parte C, la profundidad no gobierna las PC.

</details>"""

ES["ex10"] = """## Ejercicio 10 · PCA dentro de un tipo celular: respuestas a citocinas
⏱ *10 min*

Las diferencias entre tipos celulares ocultan casi todo lo demás. Un remedio habitual es analizar **un tipo celular cada vez**: seleccionamos sus muestras, elegimos sus propios genes más variables y repetimos el PCA.
Dentro de un tipo celular, lo que varía son los donantes y las citocinas.

**Qué hacer:** elige un **tipo celular** (por defecto monocitos, el más sensible) y las **PC**. Colorea por **citocina** y marca las que quieras resaltar (PBS siempre en negro), o colorea por donante o por un gen.
Los gráficos de barras de la derecha muestran los **genes con mayor peso** (loadings) en cada eje: te dicen qué *significa* el eje."""

ES["ex10_try"] = """#### ❓ Preguntas
1. En el gráfico por defecto (PC1 vs PC4 en monocitos), ¿qué citocinas se alejan de la nube principal y en qué dirección? Mira los genes detrás de cada eje.
2. Cambia los ejes a PC1 vs PC2 y colorea por **donante**. ¿Qué domina las primeras PC dentro de un tipo celular? Compara el R² de donante y citocina en la tabla.
3. Colorea por la expresión de `UTY` (un gen del cromosoma Y) y prueba varias PC: ¿qué PC separa donantes masculinos y femeninos en monocitos? (el artículo habla de 6 hombres y 6 mujeres)

<details>
<summary>💡 Respuestas</summary>

1. **IFN-β** se mueve a la derecha, a lo largo de la PC1, donde los genes principales incluyen genes inducidos por interferón (`ACOD1`, `CCL8`, `CCL7`, `IDO1`); **IL-4 e IL-13** suben a lo largo de la PC4, cuyos genes principales incluyen `ALOX15` y `WNT5B` (dianas de IL-4). **IL-4 e IL-13 caen en el mismo sitio**: tenlo presente para el Ejercicio 12.
2. **Donantes**: las diferencias basales entre personas son la mayor variación dentro de un tipo celular (R² de donante 0.72 en PC1 y 0.81 en PC2, frente a 0.18 y 0.13 de citocina). Las citocinas aparecen más claras en PC posteriores (R² 0.32 en PC3, 0.38 en PC4): los efectos individuales son menores que las diferencias entre personas. Por eso compararemos cada citocina con el PBS del **mismo donante**.
3. **PC5** (correlación 0.87 entre `UTY` y la PC5). El sexo del donante es una de las mayores fuentes de variación en monocitos, aunque nada tenga que ver con las citocinas. (En los otros tipos celulares aparece en otra PC: PC2 en T CD4, PC3 en T CD8 y NK.)

</details>"""

ES["partE"] = """---
# Parte E · Expresión diferencial: ¿qué genes cambian?

El PCA mostró que las citocinas mueven las muestras en direcciones concretas. La **expresión diferencial (DE)** hace la pregunta **gen a gen**:

> *¿La expresión del gen g con la citocina es distinta de su expresión con PBS, más de lo que esperaríamos por azar?*

Los resultados que vas a explorar se **calcularon de antemano** con **DESeq2** (ahora no se ejecuta nada): para cada uno de los 4 tipos celulares y cada una de las 87 citocinas, la comparación *citocina vs PBS* en los 12 donantes.

> 🧬 **Solo los genes inmunes.** Para que el análisis sea ligero y los gráficos legibles, la DE se hizo con **~500 genes del sistema inmune** (506 expresados): la unión de los conjuntos GO *cytokine-mediated signaling pathway* y *T cell activation* y los conjuntos Hallmark de *respuesta a interferón α / γ*. Cada comparación analiza ~450 de ellos.
> **¿Qué perdemos?** (1) La lista se construyó *a partir de genes que responden a citocinas*, así que la **fracción de DEG es mucho mayor** que en un análisis de todo el transcriptoma; no compares porcentajes. (2) **Los genes fuera de la lista no se pueden encontrar**, ni siquiera los clásicos (`ACTB`, `GBP1`, `ALOX15`, `NKG7`, `CCL17`… no están): las herramientas avisan cuando un gen no se analizó.

<details>
<summary>🔬 ¿Cómo se calcularon los resultados? (opcional)</summary>

* DESeq2 trabaja con los **conteos crudos** y los normaliza él mismo (una versión robusta del CPM).
* El **diseño** `~ donor + cytokine` es **pareado**: cada donante se compara consigo mismo, así que las diferencias entre personas se cancelan.
* Para cada gen ajusta un modelo **binomial negativo** (los conteos no son normales) y un **test de Wald** pregunta si el efecto de la citocina es mayor que el ruido.
* Como se analizan ~450 genes a la vez, los p-valores se **ajustan por comparaciones múltiples** (Benjamini–Hochberg, controlando la tasa de falsos descubrimientos, FDR): el **p-valor ajustado** es el que hay que usar.
* Los genes con muy pocos conteos (menos de 10 conteos en 12 de las 24 muestras) no se analizan.

</details>

Dos números describen cada gen:

| Número | Significado |
|---|---|
| **log₂ fold change** | el **tamaño** del cambio: +1 = se duplica (×2), +2 = ×4, −1 = se reduce a la mitad |
| **p-valor ajustado** | la **evidencia**: cuán sorprendente sería el cambio si la citocina no hiciera nada (pequeño = convincente) |

Un gen se considera **diferencialmente expresado (DEG)** cuando supera **ambos** umbrales.
"""

ES["ex11"] = """## Ejercicio 11 · Volcano plots: los genes que responden
⏱ *12 min*

Un **volcano plot** muestra todos los genes a la vez: *x* = log₂ fold change (tamaño del efecto), *y* = −log₁₀ p-valor ajustado (evidencia; más alto = más convincente).
**Rojo** = DEG sobre-expresados, **azul** = DEG infra-expresados, **gris** = no significativos. Los genes arriba y lejos del centro son los más robustos y fuertes.

**Qué hacer:** por defecto, **IFN-γ en monocitos**, con los genes `CXCL9`, `CXCL10`, `IDO1`, `STAT1` e `IRF1` resaltados. Cambia el **tipo celular** y la **citocina**, los **umbrales** y los genes a resaltar
(escribe un símbolo o elige una lista predefinida; solo se pueden encontrar los ~500 genes inmunes); las tablas muestran su fold change, p-valor ajustado y expresión media."""

ES["ex11_try"] = """#### ❓ Preguntas
1. ¿Cuántos genes suben y bajan con IFN-γ en monocitos? Mira los cinco genes resaltados: ¿cuáles son DEG y cuántas veces más altos que en PBS? ¿Alguno es significativo pero *no* es DEG?
2. Pon el umbral de fold change a 0, 1 y 2 (p-valor ajustado 0.05), y después cambia el umbral del p-valor ajustado de 0.2 a 0.001 (fold change 0.5). ¿Qué umbral importa más? ¿Qué aporta el umbral de fold change?
3. Cambia la citocina a **IL-4** y resalta `CCL26`, `CD1C`, `STAT1`, `IDO1`. ¿Qué pasa con los genes de interferón?
4. Elige **células NK** e **IL-2**. ¿Qué tienen en común los genes más sobre-expresados `IL2RA`, `CISH` y `SOCS2`?

<details>
<summary>💡 Respuestas</summary>

1. **96 suben y 95 bajan** (de 452 genes analizados). `CXCL9` (×10.5), `IDO1` (×5.5) y `CXCL10` (×3.7) suben mucho: la clásica **respuesta a IFN-γ** (reclutamiento de células T y catabolismo del triptófano). `STAT1` **no** es DEG (×1.1, p ajustado = 0.34), aunque es el factor de transcripción de esa respuesta: ya es abundante en monocitos. `IRF1` es **significativo** (p ajustado ≈ 10⁻⁵) pero su cambio es pequeño (×1.4, log₂FC 0.49 < 0.5): significativo no es lo mismo que grande.
2. Umbral de fold change 0 / 1 / 2 → **311 / 70 / 14** DEG; p ajustado 0.2 / 0.001 → 197 / 161 DEG (con fold change 0.5; 191 con 0.05). **El umbral de fold change importa mucho más**: elimina cambios estadísticamente significativos pero biológicamente pequeños, y deja los que son fiables y grandes.
3. IL-4 sobre-expresa con fuerza sus dianas (`CCL26` ×207, `CD1C` ×21) y **reduce** los genes de interferón (`IDO1` ×0.1, `STAT1` ×0.4). IL-4 e IFN-γ empujan a los monocitos en sentidos opuestos.
4. Son genes de la **señalización de IL-2**: `IL2RA` es la cadena α del receptor de alta afinidad de IL-2, y `CISH` y `SOCS2` son dianas de STAT5 que retroalimentan para apagar la señal. Las células NK responden a IL-2 a través de STAT5.

</details>
"""

ES["ex12"] = """## Ejercicio 12 · Comparar dos respuestas: ¿parecidas o distintas?
⏱ *10 min*

¿Activan dos citocinas los mismos genes? Para cada gen tenemos un log₂ fold change en la **comparación A** y en la **comparación B** (cada una = una citocina vs PBS en un tipo celular). En un diagrama de dispersión de A frente a B,
dos respuestas parecidas caen a lo largo de la **diagonal**, y el **coeficiente de correlación** lo resume: **r = 1** misma respuesta, **r = 0** sin relación, **r = −1** opuestas.
Los colores indican dónde es significativo el gen: **azul** solo en A, **naranja** solo en B, **violeta** en ambas, gris en ninguna. El panel derecho ordena las demás citocinas según lo que se parecen a A (en el mismo tipo celular).

**Qué hacer:** por defecto se comparan **IL-4 e IL-13 en monocitos**. Cambia la citocina y el tipo celular de A y B (¡pueden ser de tipos celulares distintos!), carga uno de los **ejemplos**, elige qué genes se usan para la correlación y **resalta genes**.

> ⚠️ Ambas comparaciones usan el **mismo control PBS**: el ruido de las muestras PBS acaba en las dos y puede crear o inflar una correlación. Cuidado con valores pequeños de r (≈ 0.2 o menos)."""

ES["ex12_try"] = """#### ❓ Preguntas
1. ¿Cuánto se parecen las respuestas a IL-4 e IL-13 en monocitos? Ahora cambia ambas comparaciones a **células T CD4**. ¿Qué pasa y por qué? (recuerda `IL13RA1` del Ejercicio 7)
2. Carga los ejemplos **IFN-β vs IFN-γ en monocitos** e **IFN-γ vs IL-4 en monocitos**. Compara las correlaciones y usa el panel de la derecha: ¿qué citocinas responden más como IFN-γ?
3. Carga **TNF-α en monocitos vs en células T CD4** e **IL-6 en monocitos vs en células T CD4**: ¿qué estímulo desencadena una respuesta más parecida en los dos tipos celulares? ¿Y **TNF-α vs IL-6** en monocitos?
4. Cambia *genes usados para la correlación* a "genes significativos en A o en B". ¿Cuándo cambia mucho r y cuándo no?

<details>
<summary>💡 Respuestas</summary>

1. En monocitos las dos respuestas son **casi idénticas (r = 0.97)** e IL-13 es la citocina más parecida a IL-4 de 86. En células T CD4 **r cae a 0.41** (IL-13 es solo la 3.ª). IL-4 e IL-13 comparten receptor (IL-4Rα/IL-13Rα1), pero las células T apenas expresan `IL13RA1` (~2 CPM frente a ~256 en monocitos), así que responden a IL-4 y no a IL-13. **El mismo par de estímulos puede parecer idéntico o muy distinto según el tipo celular**.
2. **IFN-β vs IFN-γ: r = 0.78** (ambos interferones activan STAT1 y los genes estimulados por interferón); **IFN-γ vs IL-4: r = −0.35**: programas opuestos (lo que uno sube, el otro lo baja). El ranking de IFN-γ en monocitos empieza con IL-12 (0.87), IFN-β e IFN-ω (0.78) e IL-15 (0.74) y termina con IL-17E y TNF-α (−0.46, −0.45).
3. **TNF-α** da una respuesta parecida en los dos tipos celulares (**r = 0.70**), **IL-6** una más específica de tipo celular (**r = 0.47**). **TNF-α vs IL-6 en monocitos: r = 0.22**, una relación débil. Mira siempre los datos antes de suponer que dos citocinas comparten un programa.
4. En monocitos, donde ambas citocinas cambian muchos genes, r apenas se mueve (IL-4 vs IL-13: 0.97 → 0.97). Donde las respuestas son más débiles, restringirse a genes que responden (células T CD4: 0.41 → 0.51; TNF-α vs IL-6 en monocitos: 0.22 → 0.28) elimina ruido y r **aumenta**. Pero también *selecciona* genes que cambian, así que la comparación deja de estar libre de sesgo: indica qué genes se usaron. Recuerda además que todo el análisis usa ya una lista de genes inmunes, enriquecida en genes que responden a citocinas: las correlaciones son mayores que en todo el genoma.

</details>
"""

ES["wrapup"] = """---
# 🏁 Cierre

Has pasado de una sola variable a comparar programas transcripcionales:

| Ejercicio | La idea a recordar |
|:--:|---|
| 1–2 | Una variable es un nombre para un valor; cada tipo de gráfico responde a una pregunta. |
| 3 | AnnData mantiene **alineados** los números (`X`), las muestras (`obs`) y los genes (`var`); `adata[muestras, genes]` selecciona un subconjunto. |
| 4 | Los conteos crudos mezclan biología con la **cantidad de material** (tamaño de librería). |
| 5 | El **CPM** divide por el tamaño de librería para que las muestras sean comparables. |
| 6 | El **logaritmo** comprime la escala, convierte cocientes en diferencias y evita que unos pocos genes dominen. |
| 7 | Los **mapas de calor** muestran patrones; solo tienen sentido *tras* normalizar. |
| 8–10 | El **PCA** es un giro que encuentra las direcciones de mayor varianza. Colorear por metadatos dice qué significa cada PC: tipo celular ≫ donante ≫ citocina. Mira **dentro de un tipo celular** para ver los efectos de las citocinas. |
| 11 | Un **DEG** necesita evidencia (**p-valor ajustado**) *y* un efecto relevante (**fold change**). El diseño pareado compara cada donante consigo mismo. |
| 12 | La **correlación de fold changes** dice si dos estímulos activan programas parecidos, y **depende del tipo celular**. |

### Para seguir
* **Enriquecimiento funcional** de las listas de DEG (¿qué rutas?): `gseapy`, `decoupler`.
* **Análisis de célula única** de los datos originales: `scanpy`; buenas prácticas en Heumos *et al.* (2023).
* El propio **Human Cytokine Dictionary** y su paquete `huCIRA` (Oesinghaus *et al.*, 2025) para interpretar tus datos con firmas de respuesta a citocinas.

### 📝 Cuestionario rápido
Tu profesor mostrará el **código QR** del cuestionario (3 preguntas, 1 minuto).

### Referencias
* Oesinghaus L. *et al.* (2025) A single-cell cytokine dictionary of human peripheral blood. *bioRxiv*, doi:10.64898/2025.12.12.693897
* Virshup I. *et al.* (2024) anndata: Access and store annotated data matrices. *J. Open Source Softw.* 9:4371
* Heumos L. *et al.* (2023) Best practices for single-cell analysis across modalities. *Nat. Rev. Genet.* 24:550–572
* Love M.I., Huber W., Anders S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.* 15:550
* Muzellec B. *et al.* (2023) PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* 39:btad547
* Benjamini Y., Hochberg Y. (1995) Controlling the false discovery rate. *J. R. Stat. Soc. B* 57:289–300
"""
