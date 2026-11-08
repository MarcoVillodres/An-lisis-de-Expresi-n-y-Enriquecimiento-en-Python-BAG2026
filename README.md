# Análisis de Expresión Diferencial y Enriquecimiento

Nombre del Estudiante: Marco Antonio Sánchez Villodres

Lenguaje elegido: Python (3.12.16)

**Repositorio:** https://github.com/<tu-usuario>/<tu-repositorio>

## 1. Instrucciones del proyecto

### Requisitos previos

- Sistema operativo **Ubuntu**.
- **Miniconda3**.
- **Conexión a internet**, gseapy consulta el servidor de Enrichr.
- Datos de entrada en `data/raw/`: `counts.csv.gz` y `metadata.csv.gz`.

### Creación del entorno

Abrir una terminal en la carpeta raíz del proyecto y ejecutar (solo la primera vez):

```bash
# Crea el entorno en ./.venv a partir de environment.yml
conda env create -p ./.venv -f environment.yml
```

El entorno instala Python 3.12, pandas, numpy, matplotlib, pydeseq2 y gseapy
con las versiones fijadas en `environment.yml`.

### Ejecución del análisis

```bash
# Dar permiso de ejecución (solo la primera vez)
chmod +x main.sh

# Ejecutar el análisis completo
./main.sh
```

`main.sh` activa el entorno local `./.venv` y ejecuta `python src/01_run_deseq2.py`. No hace falta activar el entorno a mano. El script muestra el progreso por pasos (`[1/8]` … `[8/8]`) y termina sin errores si todo va bien.

## 2. Resultados

| Ubicación | Contenido |
|---|---|
| `results/deseq2_all_results.tsv` | Resultados de todos los genes (log2FC, p-valor y padj) |
| `results/*.tsv` | Subconjuntos de genes (significativos, DE, up, down) y términos enriquecidos de GO y Reactome |
| `results/*.png` | Dispersión, volcano, MA-plot, barras top 10 y dotplots de enriquecimiento |
| `results/gen_top_padj.txt` | Gen más significativo |
| `logs/sessionInfo.txt` | Versiones de Python y de las librerías utilizadas |
| `logs/run_summary.txt` | Resumen numérico de la ejecución (comparable con el de R) |

## 3. Decisiones e incidencias

- **Enriquecimiento:** `gseapy.enrichr` trabaja con símbolos de gen (sin conversión a
  ENTREZID), con bibliotecas de Enrichr (<GO_Biological_Process_AAAA>, <Reactome_AAAA>) y
  los genes de `de_df` como universo.
- **Diferencias con R:** los resultados de enriquecimiento no coinciden exactamente
  (distinta versión y estructura de las anotaciones de GO, así como distintos filtros:
  R usa `pvalueCutoff` y `qvalueCutoff`, mientras que aquí se filtra por
  `Adjusted P-value < alpha`). 