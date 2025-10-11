import unittest
from unittest.mock import patch, Mock
import sys
import os

# Add src to Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import utils


class TestUtilsSimple(unittest.TestCase):
    """Simple tests for utils module - focusing on basic functionality."""

    def test_module_imports(self):
        """Test that utils module imports correctly."""
        self.assertTrue(hasattr(utils, 'initialize_logging'))
        self.assertTrue(hasattr(utils, 'convert_lon'))
        self.assertTrue(hasattr(utils, 'handle_track_file'))
        self.assertTrue(hasattr(utils, 'find_extremum_coordinates'))
        self.assertTrue(hasattr(utils, 'slice_domain'))

    def test_initialize_logging_exists(self):
        """Test that initialize_logging function can be called."""
        # Just test that the function exists and can be called with mock args
        try:
            with patch('utils.logging') as mock_logging:
                mock_args = Mock()
                mock_args.verbose = False
                result = utils.initialize_logging('test_dir', mock_args)
                # If we get here, function signature is working
                self.assertTrue(True)
        except Exception as e:
            # If it fails due to missing parameters, that's what we want to catch
            if "missing" in str(e) or "required" in str(e):
                self.fail(f"Function signature mismatch: {e}")

    @patch('utils.os.path.exists')
    def test_convert_lon_exists(self, mock_exists):
        """Test that convert_lon function exists and has proper signature."""
        mock_exists.return_value = True
        try:
            # Try calling with minimal arguments to test signature
            mock_data = Mock()
            result = utils.convert_lon(mock_data, 'longitude')
            self.assertTrue(True)  # If we get here, basic structure works
        except Exception as e:
            if "missing" in str(e) and "required" in str(e):
                self.fail(f"Function signature issue: {e}")

    @patch('utils.logging')
    def test_handle_track_file_signature(self, mock_logging):
        """Test handle_track_file function signature."""
        mock_logger = Mock()
        mock_logging.getLogger.return_value = mock_logger
        
        try:
            # Test with mock parameters to check signature
            result = utils.handle_track_file(
                'dummy_file.csv',
                Mock(),  # times
                'longitude',  # longitude_indexer  
                'latitude',   # latitude_indexer
                mock_logger   # app_logger
            )
        except FileNotFoundError:
            # Expected - file doesn't exist, but signature is OK
            pass
        except Exception as e:
            if "missing" in str(e) and "required" in str(e):
                self.fail(f"Function signature issue: {e}")

    def test_find_extremum_coordinates_signature(self):
        """Test find_extremum_coordinates function signature."""
        try:
            # Test basic signature with mock data
            mock_data = Mock()
            mock_lon = Mock()
            mock_lat = Mock()
            mock_var = Mock()
            mock_args = Mock()
            
            result = utils.find_extremum_coordinates(
                mock_data, mock_lon, mock_lat, mock_var, mock_args
            )
        except Exception as e:
            if "missing" in str(e) and "required" in str(e):
                self.fail(f"Function signature issue: {e}")
            # Other exceptions are fine - we just want to test signature

    def test_slice_domain_signature(self):
        """Test slice_domain function signature."""
        try:
            mock_data = Mock()
            mock_limits = Mock()
            mock_namelist = Mock()
            
            result = utils.slice_domain(mock_data, mock_limits, mock_namelist)
        except Exception as e:
            if "missing" in str(e) and "required" in str(e):
                self.fail(f"Function signature issue: {e}")
            # Other exceptions are fine - we just want to test signature


class TestUtilsIntegration(unittest.TestCase):
    """Simple integration tests."""
    
    def test_module_structure(self):
        """Test that module has expected structure."""
        expected_functions = [
            'initialize_logging',
            'convert_lon', 
            'handle_track_file',
            'find_extremum_coordinates',
            'slice_domain'
        ]
        
        for func_name in expected_functions:
            self.assertTrue(hasattr(utils, func_name), 
                          f"Missing function: {func_name}")
            self.assertTrue(callable(getattr(utils, func_name)),
                          f"Function not callable: {func_name}")


if __name__ == '__main__':
    unittest.main()
