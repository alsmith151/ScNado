import scnado
import os

def test_get_barcode_csv():
    path = scnado.get_barcode_csv()
    assert os.path.exists(path)
    assert path.endswith('barcodes_extracted.csv')

def test_data_files_exist():
    for filename in ['BC1_barcodes.txt', 'BC2_barcodes.txt', 'BC3_barcodes.txt', 'BC4_barcodes.txt']:
        path = scnado.utils.get_data_path(filename)
        assert os.path.exists(path)
