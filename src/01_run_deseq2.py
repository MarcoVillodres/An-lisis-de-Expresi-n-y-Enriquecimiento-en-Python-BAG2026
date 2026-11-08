#Definicion de los paquetes
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import gseapy as gp 
import pydeseq2 as pypy
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import os

#Creacion de las carpetas results y logs vacías
os.makedirs("results", exist_ok=True)
os.makedirs("logs", exist_ok=True)

#Definición de las variables
counts_path = "data/raw/counts.csv.gz"
metadata_path = "data/raw/metadata.csv.gz"
condition_column = "site"
plate_column = "plate"
condition_test = "Primary_Solid_Tumor"
condition_ref = "Solid_Tissue_Normal"
min_samples = 20
min_count = 10
alpha = 0.05
lfc_threshold = 1


print("[1/8] Leyendo datos de entrada (TCGA-BRCA)...")
counts_df = pd.read_csv(counts_path, index_col= 0) 
metadata = pd.read_csv(metadata_path, index_col = 0) 

counts_mat = counts_df.to_numpy() 

print(counts_df.shape) ##Comprobación de las dimensiones de la matriz generada son correctas



print("[2/8] Alineando muestras entre conteos y metadata...")
metadata = metadata.loc[counts_df.columns] 

metadata[condition_column] = metadata[condition_column].astype('category')
metadata[plate_column] = metadata[plate_column].astype('category')
metadata[condition_column] = pd.Categorical(metadata[condition_column], categories = [condition_ref, condition_test], ordered = False)

print(counts_df.columns == metadata.index) #Comprobación de que las columnas de counts_df coinciden con index.metadata


print("[3/8] Creando objeto DESeq2 y filtrando genes de baja cuenta...")
counts_df_filtered = (counts_df >= min_count).sum(axis=1) >= min_samples
counts_df_filtered = counts_df.loc[counts_df_filtered, :]
counts_df_dds = counts_df_filtered.T


dds = DeseqDataSet(counts=counts_df_dds, metadata=metadata, design=f"~ {condition_column} + {plate_column}")

print(counts_df_dds.info()) ##Comprobación de que cambian los parámetros de counts_df cambian despues del filtro, de 20502 a 17606 genes


print("[4/8] Ejecutando resultados y guardando tablas...")
dds.deseq2()
plt.ion() 
dds.plot_dispersions(save_path="results/plot_dispersion.png")
plt.close("all")


print("[5/8] Extrayendo resultados y guardando tablas...")
##Calculo de contrastes
results_crude = DeseqStats(dds, contrast=[condition_column, condition_test, condition_ref], alpha=alpha)
results_crude.summary()
results = results_crude.results_df


##Adición columna genes y ordenamiento por valores padj
print(results.head())
results["gene_id"] = results.index
print(results.head()) ##Comprobación de que aparecen la columna de los genes
results = results.sort_values("padj", na_position = "last")
print(results.head()) ##Comrpobación de que aparece ordenado según padj


##Escritura de las tablas
results.to_csv("results/deseq2_all_results.tsv", sep="\t", index=True )

##Deficion de subconjuntos de datos relevantes para el analisis de la comparación
sig_df = results[results["padj"].notna() & (results["padj"] < alpha)]
de_df = results[results["padj"].notna() & (results["padj"] < alpha) & (results["log2FoldChange"].abs() > lfc_threshold)]


##Dentro de los genes diferencialmente expresados se separan dos grupos
up_df = de_df[de_df["log2FoldChange"] > 0]
down_df = de_df[de_df["log2FoldChange"] < 0]

##Definición de los genes más significativos
top10_up_tsv = up_df.head(10)
top10_down_tsv = down_df.head(10)

##Construcción de una nueva columna definitoria de si es 10Up o 10Down
de_df["top10_label"] = "other"
de_df.loc[de_df["gene_id"].isin(top10_up_tsv["gene_id"]), "top10_label"] = "top10_up"
de_df.loc[de_df["gene_id"].isin(top10_down_tsv["gene_id"]), "top10_label"] = "top10_down"

##Generacion de las tablas conteniendo los subconjuntos de genes de interes analitico
sig_df.to_csv( "results/deseq2_sig_padj_0.05.tsv", sep = "\t", index = True)
de_df.to_csv( "results/deseq2_de_padj_0.05_lfc_1.tsv", sep = "\t", index = True)
up_df.to_csv( "results/Up_genes.tsv", sep = "\t", index = True)
down_df.to_csv( "results/Down_genes.tsv", sep = "\t", index = True)
top10_up_tsv.to_csv ( "results/top10_up_genes.tsv", sep = "\t", index = True)
top10_down_tsv.to_csv ( "results/top10_down_genes.tsv", sep = "\t", index = True)



print("[6/8] Generación de Gráficos del Análisis")

##Se generan conjuntos de datos ordenados por valores de Log2FoldChange para la generación de la Figura Top10
up_df_log2 = de_df[de_df["log2FoldChange"] > 0].sort_values("log2FoldChange", ascending=False)
down_df_log2 = de_df[de_df["log2FoldChange"] < 0].sort_values("log2FoldChange", ascending=True)
print(up_df_log2, down_df_log2) ##Comprobación de que están ordenados

##Clasificacion de los 10 genes
top10_up_log2FoldChange = up_df_log2.head(10)
top10_down_log2FoldChange = down_df_log2.head(10)


##Gráficos horizontales de los Genes top10
plt.ion()
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
axs[0].barh( y = top10_down_log2FoldChange["gene_id"], color= "#53A3FFFF", width=  top10_down_log2FoldChange["log2FoldChange"])
axs[0].set_title("Down (Tumor < Normal)")
axs[0].set_xlabel("log2FoldChange")
axs[1].barh( y = top10_up_log2FoldChange["gene_id"], color= "#FC694CFF", width=  top10_up_log2FoldChange["log2FoldChange"])
axs[1].set_title("Up (Tumor > Normal)")
axs[1].set_xlabel("log2FoldChange")
plt.tight_layout()
plt.suptitle("Top 10 Down- and Up-regulated Genes: Primary Tumor vs. Normal Tissue",
             fontsize=12, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig("results/plot_top10_genes.png", dpi=300)
plt.close("all")



##Generación del Volcano Plot
##Preparación de Datos
datos_volcano_df = results[results["padj"].notna()]

##Cribado de genes en el subset según el padj
datos_volcano_df["status"] = "No Significativo"
datos_volcano_df.loc[datos_volcano_df["gene_id"].isin(up_df["gene_id"]), "status"] = "Up"
datos_volcano_df.loc[datos_volcano_df["gene_id"].isin(down_df["gene_id"]), "status"] = "Down"

##Generación del VolcanoPlot
##Primero generacion de librería de colores segun condiciones
colors = {"Up": "#F73A60", "Down": "#004696", "No Significativo": "#76607A"}

# DataFrame para reunir los datos  y poder realizar el grafico sin bucle for
volcano_plot_df = datos_volcano_df.assign(Nlog10_padj = -np.log10(datos_volcano_df["padj"]), color = datos_volcano_df["status"].map(colors))      
plt.scatter(x = volcano_plot_df["log2FoldChange"], y = volcano_plot_df["Nlog10_padj"], s=4, c= volcano_plot_df["color"])

##Líneas de los umbrales según valores de lfc_threshold
plt.axvline(lfc_threshold, linestyle="--", color = "black", linewidth= 0.8)
plt.axvline(-lfc_threshold, linestyle = "--", color ="black", linewidth = 0.8)
plt.axhline(-np.log10(alpha), linestyle="--", color="black", linewidth=0.8)

##Títulos
plt.title("Volcano plot: Primary Solid Tumor vs. Solid Tissue Normal")
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10(padj)")


##Realizo subset del dataframe preparado de volcano plot con solo los genes mas significativos para nombrarlos en el grafico
top10_label_forVolcano_df = volcano_plot_df.head(10)

##Adicion de los nombres de los 10 genes mas significativos, .....  con bucle for
for _, row in top10_label_forVolcano_df.iterrows():
    plt.text(row["log2FoldChange"], row["Nlog10_padj"] + 1 , row["gene_id"], ha="center", fontsize=5)

##Adición de rejillas
plt.grid(True, linewidth=0.3, color="lightgrey")
plt.savefig("results/plot_volcano.png", dpi=300, bbox_inches="tight")
plt.close()


##Generación del plot MA

##No se puede realizar sobre el conjunto de datos resultante de DatSetStats como en R puestoq ue ni es un cunjunto en dataframe, 
##ni esta con los genes como nombres de las filas, ni ordenado por valores de padj, ni sin los valoresNA de padj
## y por ultimo no contiene la columna status desde la cual se determina el codigo de colores, se realizará de 
##los datos ordenados de volcano_df, se conservará el código de colores del volcano plot.

plt.scatter(x = np.log10(datos_volcano_df["baseMean"]), y = datos_volcano_df["log2FoldChange"], s = 4, c = datos_volcano_df["status"].map(colors))

##Adición del eje centrado en 0
plt.axhline(0, linestyle="--", color="black", linewidth=0.8)

plt.title("DESeq2 MA plot")
plt.xlabel("Mean of normalized counts log10(baseMean)")
plt.ylabel("log2FoldChange")

##Adición de la rejilla
plt.grid(True, linewidth=0.3, color="lightgrey")

plt.savefig("results/plot_ma.png", dpi=300, bbox_inches="tight")
plt.close()


##Se muestra el gen más signifitcativo, se incluyó el paso 5.7 junto con el 5.6 en el paso [6/8]
gen_top_padj = f"Gen cuya variación de la expresión es más significativo: {de_df.iloc[0]['gene_id']} que presenta un log2FoldChange = {de_df.iloc[0]['log2FoldChange']} y un padj = {de_df.iloc[0]['padj']}"
open("results/top_ge.txt", "w").write(gen_top_padj)


print("[7/8] Enriquecimiento GO y Reactome...")

##Definición de los conjuntos de datos a usar en el Enriquecimiento desde el conjunto de datos de de_df correspondiente en python
universo = de_df["gene_id"].unique().tolist()
genes_Up = up_df["gene_id"].unique().tolist()
genes_Down = down_df["gene_id"].unique().tolist()


##Antes de realizar el análisis GO y Reactome vemos que paquetes son los más adecuados a utilizar dentro de
##los de humano para realizar el análsis
print([ñ for ñ in gp.get_library_name(organism="human") if "GO_Biological" in ñ or "Reactome" in ñ])

##Se observan que las librerias de genes de GO biological process y Reactome más recientes
##datan de 2026 y de 2024 respectivamente, al estar más actualizadas se escogerán susodichas para el análisis.

##Tras buscar información en la base de datos de gseapy se observa que el método de gp.enrich utiliza el método de Benjamini-Hochberg por defecto 
##y que no presenta la posibilidad filtrar los datos por q-valores.

##Análisis de GO en genes Up y Down mediante función
def gseapy_enrich(Name_title, Name_file, genes_constrast, gene_Set, universo = universo, alpha = alpha):
    enricHH = gp.enrichr(gene_list = genes_constrast, gene_sets = gene_Set, organism = "human", background = universo, outdir = None)

    ##Cribado de los datos otenidos según padj
    enricHH_results =  enricHH.results[enricHH.results["Adjusted P-value"] < alpha]

    ##TSV y PNG
    enricHH_results.to_csv(f"results/{Name_file}.tsv", sep="\t", index=False)
    gp.dotplot(enricHH_results, column = "Adjusted P-value", top_term = 10, title = f"{Name_title} Significative", ofname=f"results/plot_{Name_file}.png")  
    return(enricHH_results)


##Llamadas de la función
gseapy_enrich("Go BP Up1","go_bp_up", genes_Up, "GO_Biological_Process_2026")
gseapy_enrich("Go BP Down1","go_bp_down", genes_Down, "GO_Biological_Process_2026")
gseapy_enrich("Reactome Up1","reactome_up", genes_Up, "Reactome_Pathways_2024")
gseapy_enrich("Reactome Down1","reactome_down", genes_Down, "Reactome_Pathways_2024")


print("[8/8] Guardando la trazabilidad de la sesión.....")

##Definición de paquetes utilizados y sus versiones
import sys
definitoria_paquetes = [ f"Python: {sys.version}", f"numpy: {np.__version__}", f"pandas: {pd.__version__}",
                         f"matplotlib: {plt.matplotlib.__version__}", f"os: {sys.version}",
                         f"gseapy: {gp.__version__}", f"pydeseq.dds: {pypy.__version__}", f"sys: {sys.version}"]
open("logs/sessionInfo.txt", "w").write("\n".join(definitoria_paquetes))

##Definición de los parámetros que conforman el summary
summary_lines = [f"n_samples: {counts_df.shape[1]}",f"n_genes_input: {counts_df.shape[0]}",
                 f"n_genes_after_filter: {dds.n_vars}",f"n_significant_padj_lt_{alpha}: {len(sig_df)}",
                 f"n_de_absLog2FC_gt_{lfc_threshold}_padj_lt{alpha}: {len(de_df)}",f"n_up: {len(up_df)}",
                 f"n_down: {len(down_df)}"]
open("logs/run_summary.txt", "w").write("\n".join(summary_lines))


