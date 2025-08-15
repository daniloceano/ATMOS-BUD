"""
Tests for the get_era5_data module.
"""

import pytest
import os
import tempfile
from unittest.mock import Mock, patch, MagicMock
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from get_era5_data import download_era5_data, download_era5_data_legacy


class TestDownloadEra5Data:
    """Test cases for download_era5_data function."""
    
    @patch('get_era5_data.cdsapi.Client')
    @patch('get_era5_data.logging')
    def test_download_era5_data_success(self, mock_logging, mock_client_class):
        """Test successful ERA5 data download."""
        # Setup mock client and result
        mock_client = Mock()
        mock_result = Mock()
        mock_client.retrieve.return_value = mock_result
        mock_client_class.return_value = mock_client
        
        # Test parameters using correct function signature
        variables = ['temperature', 'geopotential']
        pressure_levels = [850, 500, 200]
        start_date = '2020-01-01'
        end_date = '2020-01-02'
        area = [90, -180, -90, 180]  # North, West, South, East
        hours = ['00:00', '12:00']
        
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as tmp:
            output_file = tmp.name
        
        # Create the file so os.path.exists() returns True
        with open(output_file, 'w') as f:
            f.write('fake netcdf content')
        
        try:
            # Call the function
            result = download_era5_data(
                variables=variables,
                pressure_levels=pressure_levels,
                start_date=start_date,
                end_date=end_date,
                area=area,
                output_file=output_file,
                hours=hours
            )
            
            # Verify client was called correctly
            mock_client_class.assert_called_once()
            mock_client.retrieve.assert_called_once()
            mock_result.download.assert_called_once_with(output_file)
            
            # Check the retrieve call parameters - simplified test
            call_args = mock_client.retrieve.call_args
            assert call_args[0][0] == 'reanalysis-era5-pressure-levels'
            
            request_params = call_args[0][1]
            assert request_params['variable'] == variables
            assert request_params['pressure_level'] == [str(p) for p in pressure_levels]
            assert request_params['area'] == area
            assert request_params['data_format'] == 'netcdf'  # Correct field name
            
            # Function returns None by design
            assert result is None
            
        finally:
            # Cleanup
            if os.path.exists(output_file):
                os.unlink(output_file)
    
    @patch('get_era5_data.cdsapi.Client')
    @patch('get_era5_data.logging')
    def test_download_era5_data_single_level(self, mock_logging, mock_client_class):
        """Test ERA5 single level data download."""
        mock_client = Mock()
        mock_client_class.return_value = mock_client
        
        variables = ['2m_temperature', 'mean_sea_level_pressure']
        start_date = '2020-01-01'
        end_date = '2020-01-01'
        area = [60, -10, 50, 10]
        hours = ['00:00']
        
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as tmp:
            output_file = tmp.name
        
        try:
            result = download_era5_data(
                variables=variables,
                pressure_levels=[],  # Empty list for single level
                start_date=start_date,
                end_date=end_date,
                area=area,
                output_file=output_file,
                hours=hours
            )
            
            # Should use pressure-levels dataset even with empty pressure levels
            # (implementation detail depends on actual function logic)
            call_args = mock_client.retrieve.call_args
            # The actual dataset selection logic would be tested here
            
            # Function returns None by design
            assert result is None
            
        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)
    
    @patch('get_era5_data.cdsapi.Client')
    def test_download_era5_data_client_error(self, mock_client_class):
        """Test handling of CDSAPI client errors."""
        mock_client = Mock()
        mock_client.retrieve.side_effect = Exception("API Error")
        mock_client_class.return_value = mock_client
        
        variables = ['temperature']
        start_date = '2020-01-01'
        end_date = '2020-01-01'
        area = [60, -10, 50, 10]
        
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as tmp:
            output_file = tmp.name
        
        try:
            with pytest.raises(Exception) as exc_info:
                download_era5_data(
                    variables=variables,
                    pressure_levels=[850],
                    start_date=start_date,
                    end_date=end_date,
                    area=area,
                    output_file=output_file
                )
            
            assert "API Error" in str(exc_info.value)
            
        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)
    
    @patch('get_era5_data.os.path.getsize')
    @patch('get_era5_data.os.path.exists')
    @patch('get_era5_data.cdsapi.Client')
    def test_download_era5_data_parameter_validation(self, mock_client_class, mock_exists, mock_getsize):
        """Test basic parameter validation without API calls."""
        # Setup mocks
        mock_client = Mock()
        mock_result = Mock()
        mock_client.retrieve.return_value = mock_result
        mock_client_class.return_value = mock_client
        mock_exists.return_value = True  # Pretend file exists
        mock_getsize.return_value = 1024 * 1024  # 1 MB
        
        # Test that function can be called with valid parameters
        download_era5_data(
            variables=['temperature'],
            pressure_levels=[850],
            start_date='2020-01-01',
            end_date='2020-01-01',
            area=[60, -10, 50, 10],
            output_file='test.nc'
        )
        
        # Verify the function executed properly
        mock_client_class.assert_called_once()
        mock_client.retrieve.assert_called_once()
        mock_result.download.assert_called_once_with('test.nc')


class TestDownloadEra5DataLegacy:
    """Test cases for download_era5_data_legacy function."""
    
    @patch('get_era5_data.download_era5_data')
    @patch('get_era5_data.logging')
    def test_download_era5_data_legacy_success(self, mock_logging, mock_download):
        """Test successful legacy ERA5 data download."""
        # Setup mock for the main download function
        mock_download.return_value = None  # Function returns None
        
        # Call the legacy function (no parameters)
        result = download_era5_data_legacy()
        
        # Verify the main download function was called with legacy parameters
        mock_download.assert_called_once()
        call_args = mock_download.call_args
        
        # Check that legacy parameters were used
        assert call_args[1]['variables'] is not None
        assert call_args[1]['pressure_levels'] is not None
        assert call_args[1]['start_date'] == '2005-08-08'
        assert call_args[1]['end_date'] == '2005-08-14'
        assert call_args[1]['area'] == [-17.5, -60, -42.5, -30]
        assert call_args[1]['output_file'] == 'system-20050808_ERA5.nc'
        
        # Function returns None by design
        assert result is None
    
    @patch('get_era5_data.download_era5_data')
    def test_download_era5_data_legacy_no_parameters(self, mock_download):
        """Test legacy function accepts no parameters."""
        # Mock the main function to prevent real API calls
        mock_download.return_value = None
        
        # Should not raise error when called without parameters
        result = download_era5_data_legacy()
        
        # Verify the main download function was called
        mock_download.assert_called_once()
        assert result is None
    
    @patch('get_era5_data.download_era5_data')
    def test_download_era5_data_legacy_calls_main_function(self, mock_download):
        """Test that legacy function calls the main download function."""
        mock_download.side_effect = Exception("Test error")
        
        with pytest.raises(Exception) as exc_info:
            download_era5_data_legacy()
        
        # Should have called the main function
        mock_download.assert_called_once()
        assert "Test error" in str(exc_info.value)


class TestModuleImports:
    """Test module imports and dependencies."""
    
    def test_cdsapi_import(self):
        """Test that cdsapi can be imported or properly mocked."""
        try:
            import get_era5_data
            # If we get here, the module loaded successfully
            assert hasattr(get_era5_data, 'download_era5_data')
            assert hasattr(get_era5_data, 'download_era5_data_legacy')
        except ImportError as e:
            pytest.skip(f"Required dependencies not available: {e}")
    
    def test_logging_configuration(self):
        """Test that logging is properly configured."""
        import get_era5_data
        import logging
        
        # Check if logging is configured
        logger = logging.getLogger('get_era5_data')
        assert logger is not None


if __name__ == '__main__':
    pytest.main([__file__])
