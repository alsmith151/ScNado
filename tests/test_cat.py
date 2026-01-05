import pytest
from unittest.mock import MagicMock, patch
import numpy as np
from scnado.cat import process_cat

def test_process_cat_mock():
    # Mock snapatac2
    with patch('snapatac2.pp.import_fragments') as mock_import, \
         patch('snapatac2.metrics.tsse') as mock_tsse, \
         patch('snapatac2.pp.add_tile_matrix') as mock_tile, \
         patch('snapatac2.pp.select_features') as mock_select:
        
        # Setup mock data
        mock_adata = MagicMock()
        mock_adata.X = np.array([[1, 0], [0, 1]])
        mock_adata.__len__.return_value = 1
        mock_adata.__getitem__.return_value = mock_adata
        
        mock_import.return_value = mock_adata
        
        # Call the function
        result = process_cat(['dummy.bed'], 'output.h5ad')
        
        # Verify calls
        mock_import.assert_called_once()
        mock_tsse.assert_called_once()
        mock_tile.assert_called_once()
        mock_select.assert_called_once()
        
        assert len(result) == 1
        assert result[0] == mock_adata
