from .utils import get_barcode_csv, get_data_path

__version__ = "0.1.0"


def get_default_barcode_csv() -> str:
    """Return the path to the packaged barcode reference CSV."""
    return get_barcode_csv()


def get_barcode_whitelists() -> list[str]:
    """Return BC3, BC2, BC1, BC4 whitelist paths in the order STARsolo expects."""
    return [
        get_data_path("BC3_barcodes.txt"),
        get_data_path("BC2_barcodes.txt"),
        get_data_path("BC1_barcodes.txt"),
        get_data_path("BC4_barcodes.txt"),
    ]
