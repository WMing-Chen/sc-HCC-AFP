import anndata
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy.external as sce
sc.settings.verbosity = 3 
# sc.logging.print_versions() 
sc.settings.set_figure_params(dpi=100) 
import os 

results_file = 'data/HCC_scRNA_seq.h5ad' 
data = sc.read_h5ad("data/HCC_scRNA.h5ad")
data.var_names_make_unique() 

sc.pl.highest_expr_genes(data, n_top=20)
sc.pl.violin(data, ['nFeature_RNA', 'nCount_RNA', 'percent.mt'],jitter=0, multi_panel=True, show=False)


sc.pp.normalize_total(data, target_sum=1e4) ##
sc.pp.log1p(data)
data.raw = data
sc.pp.highly_variable_genes(data, n_top_genes=2000) #, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(data, show=False)

plt.savefig("Fig1/highly_variable_genes.pdf")
data = data[:, data.var.highly_variable] 

sc.pp.scale(data, max_value=10)
sc.tl.pca(data, svd_solver='arpack')# svd_solver   SVD 
sc.pl.pca_variance_ratio(data, log=True, show = False)
plt.tight_layout()
plt.savefig("Fig1/PCA.pdf")
sce.pp.harmony_integrate(data, 'Sample') # 
sc.pp.neighbors(data, use_rep = 'X_pca_harmony', n_pcs=50) # 
sc.tl.louvain(data, resolution=1) #sc.tl.louvain(data)  sc.tl.leiden(data) 
sc.tl.umap(data)

sc.tl.louvain(data, key_added="louvain_res1_3", resolution=1.3)
sc.tl.umap(data)
data.write(results_file)

data.obs['Cluster'] = data.obs['louvain_res1_3']

cols = ["#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
          "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
          "#33A02C", "#B2DF8A", "#55B1B1", "#8DD3C7", "#A6761D",
          "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
          "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
          "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00", "#28B461"]

fig, ax = plt.subplots(figsize=(7.5, 6))
sc.pl.umap(data, color=['Cluster'], ax=ax, show=False, size=5, palette=cols, title = 'Cluster') # , legend_loc='on data'
plt.tight_layout()
# plt.savefig("Fig1/Umap.pdf")


from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('custom', ["lightgrey", "orange", "red"])
sc.pl.umap(data, 
           color=["CD3D", "CD8A", "TRDC", # CD8+ T, γδT, 
              "IL7R", "CD3E", "CCR6",     #CD4+ T, T reg,"CD4", "FOXP3", "CD14", "CD19", "IL2RA", "ITGB1", "CD3E", 
              "NKG7", "GNLY", "FCGR3A", # NK  CD56:"NCAM1",  CD16:"FCGR3A"
              "CD79A", "MS4A1", # B
              "IGHG1", "MZB1"  # plasma cells
                   ],
           cmap=cmap,
           show=False)
sc.pl.umap(data, 
           color=[ "CD14", "FCN1", "APOBEC3A",  # monocytes, CD16:"FCGR3A"  CD11C:ITGAX  "THBS1", 
              "CD68", "CD163", # macrophages, 
              "CD1C", "LAMP3", "FCER1A", # dendritic cells, 
              "TPSAB1" # mast cells 
              # "CSF3R", "S100A8" # neutrophils
                   ], 
           cmap=cmap,
           show=False)
sc.pl.umap(data, 
           color=[
              "VWF",  # endothelial cells   PECAM1
              "COL1A1", "COL1A2", "DCN", "LUM", # fibroblasts
              "KRT8", "EPCAM", "KRT7", "KRT19" # epithelial cells
                   ], 
           cmap=cmap,
           vmax=2,
           show=False)
# create a dictionary to map cluster to annotation label
cluster2annotation = {
     '5': 'NK cells', '9': 'NK cells',
     '0': 'CD4+ T cells', '6': 'CD4+ T cells', '15': 'CD4+ T cells', '25': 'CD4+ T cells','23': 'CD4+ T cells',
     '1': 'CD8+ T cells', '11': 'CD8+ T cells', '14': 'CD8+ T cells', 
     '12': 'B cells', '18': 'Plasma cells',
     '3': 'Endothelial', '16': 'Endothelial', '26': 'Endothelial', '24': 'Endothelial', '30': 'Endothelial',
     '8': 'Mesenchymal cells', '17': 'Mesenchymal cells',
     '27':'Mast cells',
     '10': 'Monocyte', 
     '13': 'DC', '21': 'DC', 
     '2': 'Macrophage', '20': 'Macrophage', '22': 'Macrophage', '28': 'Macrophage',
     '4': 'Tumor', '29': 'Tumor', '7': 'Tumor', '19': 'Tumor' 
}
data.obs['Major_Celltype'] = data.obs['Cluster'].map(cluster2annotation).astype('category')
cols = {"CD4+ T cells": "#E78AC3", "CD8+ T cells": "#C6EE8F", "NK cells": "#28B461",
        "B cells": "#aa8282", "Plasma cells": "#d4b7b7", "Macrophage": "#55B1B1", "DC": "#36489E", "Monocyte": "#B2DCEE", 
        "Tumor": "#DA2917", "Mast cells": "#FA7E4D", 
        "Endothelial": "#A680B9", "Mesenchymal cells": "#F0C674"
        }
fig, ax = plt.subplots(figsize=(7, 5))
sc.pl.umap(data, color=['Major_Celltype'], ax=ax, show=False, palette=cols, size=5)

plt.tight_layout()
plt.savefig("Fig1/Umap_Major_Celltype.pdf")
data.write(results_file)

marker_genes = ["IL7R", "CD3E", "CCR6", #CD4+ T, T reg, "CD4", "CCR6", "TRBC2", "CCR7", "FOXP3", 
                "CD3D", "CD8A", # CD8+ T
                "NKG7", "GNLY", # NK  CD56:NCAM1   CD16:"FCGR3A",  "NCAM1", 
                "CD79A", "MS4A1", # B
                "IGHG1", "MZB1",  # plasma cell
                "CD68", "CD163", # Macrophage, 
                "CD14", "FCN1",  # Monocyte, CD16:"FCGR3A"  CD11C:ITGAX  "THBS1", , "APOBEC3A"
                "CD1C", "FCER1A", # DC,  "LAMP3", 
                "COL1A1", "COL1A2",  # Mesenchymal cells fibroblasts "DCN",
                "TPSAB1", # Mast cells
                "VWF",  # endothelial cells  
                "KRT8", "GPC3"# Tumor epithelial cells , "EPCAM", "KRT7""KRT19"
               ]
new_cluster_names = ['CD4+ T cells', 'CD8+ T cells', 'NK cells', 'B cells', 'Plasma cells',
                     'Macrophage', 'Monocyte', 'DC', 
                     'Mesenchymal cells','Mast cells', 'Endothelial',
                     'Tumor']
data.obs['Major_Celltype'] = pd.Categorical(data.obs['Major_Celltype'], categories=new_cluster_names, ordered=True)

fig, ax = plt.subplots(figsize=(10, 5)) 
sc.pl.dotplot(data, marker_genes, ax=ax, groupby='Major_Celltype', show=False)
plt.tight_layout()
plt.savefig("Fig1/gene_Dot_celltype.pdf")


