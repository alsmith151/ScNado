from . import cat
from . import rna
from . import multiome
from .utils import get_barcode_csv, get_data_path

__version__ = "0.1.0"


def get_default_barcode_csv() -> str:
    """Return the path to the packaged barcode reference CSV."""
    return get_barcode_csv()
