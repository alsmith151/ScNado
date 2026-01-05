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

def process_rna(sample_dirs, sample_names, output_h5ad, min_genes=100, min_cells=3, n_top_genes=2000):
    """
    Process RNA data from multiple samples.
    """
    adatas = []
    for d, name in zip(sample_dirs, sample_names):
        path = Path(d)
        # Handle different possible STARsolo output structures
        raw_path = path / 'GeneFull_Ex50pAS' / 'raw'
        if not raw_path.exists():
            raw_path = path / 'Gene' / 'raw'
            
        adata = read_dataset(raw_path)
        adata.obs['sample'] = name
        adatas.append(adata)
        
    if len(adatas) > 1:
        adata = sc.concat(adatas, label='sample', keys=sample_names)
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
    
    # HVG and DimRed
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key="sample")
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
    
    # Condition from barcode (specific to this experiment's barcode naming)
    # In the notebook: [bc.split('_')[-1].replace('-1', '') for bc in adata.obs_names]
    # We might want to make this more flexible or keep it as is for now.
    adata.obs['condition'] = [bc.split('_')[-1].replace('-1', '') for bc in adata.obs_names]
    
    adata.write(output_h5ad)
    return adata
