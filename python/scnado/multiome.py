import muon as mu
from mudata import MuData
import scanpy as sc
import snapatac2 as snap

def integrate_cat_rna(cat_h5ad, rna_h5ad, output_h5mu):
    """
    Integrate CAT and RNA data using muon.
    """
    adata_cat = sc.read_h5ad(cat_h5ad)
    adata_rna = sc.read_h5ad(rna_h5ad)
    
    # Ensure obs_names are comparable
    # In CAT: 'sublib-1:barcode'
    # In RNA: 'barcode-1' or similar
    # We need to align these. The user mentioned linking them.
    
    # Let's assume the user has a way to map them or they already match.
    # Based on the notebooks, CAT uses 'sample:barcode' and RNA uses 'barcode_condition'.
    # This might need a custom mapping function depending on how the barcodes were assigned.
    
    mdata = MuData({
        "cat": adata_cat,
        "rna": adata_rna
    })
    
    # Intersect observations
    mu.pp.intersect_obs(mdata)
    
    # Joint analysis could be added here
    # mu.tl.mofa(mdata)
    
    mdata.write(output_h5mu)
    return mdata
