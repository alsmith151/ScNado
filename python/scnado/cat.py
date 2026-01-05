import snapatac2 as snap
import numpy as np
from pathlib import Path
import scanpy as sc

def process_cat(fragment_files, output_h5ad, chrom_sizes=snap.genome.hg38, min_num_fragments=500, bin_size=5000, n_features=50000):
    """
    Process CUT&TAG fragments into a processed AnnData object.
    """
    # Import fragments
    data = snap.pp.import_fragments(
        fragment_files,
        chrom_sizes=chrom_sizes,
        sorted_by_barcode=False,
        min_num_fragments=min_num_fragments
    )
    
    # Basic QC and tile matrix
    snap.metrics.tsse(data, chrom_sizes)
    snap.pp.add_tile_matrix(data, bin_size=bin_size)
    snap.pp.select_features(data, n_features=n_features)
    
    # Binarize
    for i in range(len(data)):
        data[i].X = (data[i].X > 0).astype(np.float32)
    
    # If multiple files, they are returned as a list of AnnData objects by import_fragments
    # but we want to save them. If it's a single file, it might be just one AnnData.
    if isinstance(data, list):
        # For simplicity in the pipeline, we might want to handle one at a time or concat
        # The notebook concats them into an AnnDataSet later.
        # Here we'll just return the list or handle the first one if that's what's expected.
        # In Snakemake, we'll likely call this per sample.
        return data
    return [data]

def analyze_cat_dataset(h5ad_files, output_dataset, n_features=50000, resolution=0.5):
    """
    Analyze a collection of processed CAT AnnData objects.
    """
    data = snap.AnnDataSet(
        adatas=[(Path(p).stem, p) for p in h5ad_files],
        filename=output_dataset
    )
    
    # Unique cell IDs
    unique_cell_ids = [sa + ':' + bc for sa, bc in zip(data.obs['sample'], data.obs_names)]
    data.obs_names = unique_cell_ids
    
    # Condition from barcode
    data.obs['condition'] = [bc.split('-')[-1] for bc in data.obs_names]
    
    # Dimensionality reduction and clustering
    snap.pp.select_features(data, n_features=n_features)
    snap.tl.spectral(data)
    snap.pp.knn(data)
    snap.tl.leiden(data, resolution=resolution)
    snap.tl.umap(data)
    
    # Peak calling
    snap.tl.macs3(data, groupby='leiden', replicate='sample')
    
    return data

def export_coverage(data, out_dir, groupby='condition'):
    """
    Export bigWig coverage tracks.
    """
    snap.ex.export_coverage(
        data,
        groupby=groupby,
        out_dir=out_dir,
    )
