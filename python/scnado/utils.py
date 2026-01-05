import importlib.resources as pkg_resources
from . import data

def get_data_path(filename):
    """Get the path to a file in the scnado/data directory using package resources."""
    return str(pkg_resources.files(data).joinpath(filename))

def get_barcode_csv():
    """Get the path to the default barcodes_extracted.csv file."""
    return get_data_path('barcodes_extracted.csv')
