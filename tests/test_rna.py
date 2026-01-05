import pytest
import scnado.rna as rna
import pandas as pd
import numpy as np
from pathlib import Path
import scanpy as sc
import scipy.sparse as sp

def test_read_dataset(tmp_path):
    # Create a mock STARsolo output directory
    d = tmp_path / "raw"
    d.mkdir()
    
    # Create matrix.mtx
    # Format: 
    # %%MatrixMarket matrix coordinate real general
    # %
    # n_genes n_barcodes n_entries
    # gene_idx barcode_idx value
    n_genes = 10
    n_barcodes = 5
    matrix = sp.coo_matrix(np.random.randint(0, 100, (n_genes, n_barcodes)))
    
    with open(d / "matrix.mtx", "w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write("%\n")
        f.write(f"{n_genes} {n_barcodes} {matrix.nnz}\n")
        for r, c, v in zip(matrix.row, matrix.col, matrix.data):
            f.write(f"{r+1} {c+1} {v}\n")
            
    # Create barcodes.tsv
    barcodes = [f"BC-{i}" for i in range(n_barcodes)]
    pd.DataFrame(barcodes).to_csv(d / "barcodes.tsv", sep="\t", index=False, header=False)
    
    # Create features.tsv
    features = [[f"ENSG{i}", f"Gene{i}", "Gene Expression"] for i in range(n_genes)]
    pd.DataFrame(features).to_csv(d / "features.tsv", sep="\t", index=False, header=False)
    
    # Test read_dataset
    adata = rna.read_dataset(d)
    
    assert adata.shape == (n_barcodes, n_genes)
    assert list(adata.obs_names) == barcodes
    assert list(adata.var_names) == [f"Gene{i}" for i in range(n_genes)]
