This is scripts can be used for converting `scanpy` object to  `Seurat` object especially you get an error while using `ReadH5AD`

Its a work arround using `reticulate`, mostly modified from [here](https://theislab.github.io/scanpy-in-R/).
```
suppressPackageStartupMessages({
    library("reticulate")
    library("ggplot2")
    library("scater")
    library("Seurat")
})
```
## Import scanpy and object
``` 
sc <- import("scanpy")
adata <- sc$read_h5ad("myAnndata.h5ad")
adata
````
```
AnnData object with n_obs × n_vars = 21391 × 30983
    obs: 'n_genes', 'n_counts', 'percent_mito', 'leiden_res_1.3', 'leiden_res_2',  'leiden_labels', 'leiden_res_0.9', 'leiden_res_1.1', 'leiden_res_0.5'
    var: 'n_cells', 'percent_cells', 'robust', 'highly_variable_features', 'mean', 'var', 'hvf_loess', 'hvf_rank'
    uns: 'Channels', 'PCs', 'fmat_highly_variable_features', 'leiden_labels_colors', 'leiden_res_0.5_colors', 'leiden_res_0.9_colors', 'leiden_res_1.1_colors', 'pca'
    obsm: 'X_pca',  'X_umap'
    varm: 'PCs', 'de_res'
    layers: 'counts', 'winsorized'
    obsp: 'connectivities', 'distances'
```
## Format count and normalized data for their respective slots
```
counts <- t(adata$layers["counts"])
colnames(counts) <- adata$obs_names$to_list()
rownames(counts) <- adata$var_names$to_list()
counts <- Matrix::Matrix(as.matrix(counts), sparse = T)

data <- t(adata$layers["winsorized"])
colnames(data) <- adata$obs_names$to_list()
rownames(data) <- adata$var_names$to_list()
data <- Matrix::Matrix(as.matrix(data), sparse = T)
```
## Create the Seurat object
```
seurat <- CreateSeuratObject(counts)
seurat <- SetAssayData(seurat, "data", data)
seurat <- AddMetaData(seurat, adata$obs)
```
## Import dim. reductions
```
embedding <- adata$obsm["X_umap"]
rownames(embedding) <- adata$obs_names$to_list()
colnames(embedding) <- c("umap_1", "umap_2")
seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")
seurat <- SetIdent(seurat, value = "leiden_res_0.5")
seurat
```
```
An object of class Seurat 
30983 features across 21391 samples within 1 assay 
Active assay: RNA (30983 features, 0 variable features)
 1 dimensional reduction calculated: umap
```

Now you can use it for down stream analysis and visualisation
