import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

def read_dataset(path: Path) -> sc.AnnData:
    """
    Read STARsolo output matrices.
    """
    adata = sc.read_mtx(path / 'matrix.mtx').T
    barcodes = pd.read_csv(path / 'barcodes.tsv', header=None, sep='\t')
    genes = pd.read_csv(path / 'features.tsv', header=None, sep='\t')
    adata.obs_names = barcodes[0].values
    adata.var_names = genes[1].values
    return adata

def process_rna(filtered_dirs, output_h5ad, min_genes=100, min_cells=3, n_top_genes=2000):
    """
    Process RNA data from one or more STARsolo filtered matrix directories.

    Each entry in filtered_dirs should be the path to a STARsolo filtered/
    directory (i.e. Solo.out/GeneFull_Ex50pAS/filtered/).
    """
    adatas = []
    sample_names = []
    for d in filtered_dirs:
        path = Path(d)
        name = path.parent.parent.parent.name  # {sample}_{sublibrary} from Solo.out grandparent
        adata = read_dataset(path)
        adata.obs['sublibrary'] = name
        adatas.append(adata)
        sample_names.append(name)
        
    if len(adatas) > 1:
        adata = sc.concat(adatas, label='sublibrary', keys=sample_names)
    else:
        adata = adatas[0]

    adata.obs_names_make_unique()

    # QC
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    # Filtering
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Normalization
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # HVG and dimensionality reduction
    batch_key = "sublibrary" if len(adatas) > 1 else None
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

    adata.write(output_h5ad)
    return adata
