import unittest
from unittest.mock import patch
import sys
import os
import pytest

# Add src to Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from cli_interface import parse_arguments


class TestParseArgumentsFixed(unittest.TestCase):
    """Fixed tests for parse_arguments function without improper mocking."""

    def test_parse_arguments_missing_framework(self):
        """Test parse_arguments when required framework argument is missing."""
        test_args = ['input_file.nc']
        
        # Should exit due to missing required framework argument
        with pytest.raises(SystemExit):
            parse_arguments(test_args)
    
    def test_parse_arguments_multiple_frameworks(self):
        """Test parse_arguments with multiple mutually exclusive frameworks."""
        test_args = [
            'input_file.nc',
            '--choose',
            '--fixed'
        ]
        
        # This should fail because only one framework is allowed
        with pytest.raises(SystemExit):
            parse_arguments(test_args)
    
    def test_parse_arguments_default_values(self):
        """Test parse_arguments with default values."""
        test_args = ['input_file.nc', '--fixed']
        
        args = parse_arguments(test_args)
        
        # Check defaults
        self.assertEqual(args.infile, 'input_file.nc')
        self.assertTrue(args.fixed)
        self.assertFalse(args.track)
        self.assertFalse(args.choose)
        self.assertFalse(args.verbose)
        self.assertEqual(args.level, 850)  # default
        self.assertEqual(args.track_vorticity, 'min')  # default
        self.assertEqual(args.track_geopotential, 'min')  # default
        self.assertFalse(args.gfs)
        self.assertFalse(args.outname)
        self.assertTrue(args.save_nc_file)  # default True
    
    def test_parse_arguments_verbose_mode(self):
        """Test verbose mode parsing."""
        test_args = ['input_file.nc', '--verbose', '--track']
        
        args = parse_arguments(test_args)
        
        self.assertTrue(args.verbose)
        self.assertTrue(args.track)
        self.assertFalse(args.fixed)
        self.assertFalse(args.choose)
    
    def test_parse_arguments_framework_selection(self):
        """Test framework selection arguments."""
        # Test choose framework
        test_args_choose = ['input_file.nc', '--choose']
        args_choose = parse_arguments(test_args_choose)
        self.assertTrue(args_choose.choose)
        self.assertFalse(args_choose.fixed)
        self.assertFalse(args_choose.track)
        
        # Test fixed framework
        test_args_fixed = ['input_file.nc', '--fixed']
        args_fixed = parse_arguments(test_args_fixed)
        self.assertTrue(args_fixed.fixed)
        self.assertFalse(args_fixed.choose)
        self.assertFalse(args_fixed.track)
        
        # Test track framework
        test_args_track = ['input_file.nc', '--track']
        args_track = parse_arguments(test_args_track)
        self.assertTrue(args_track.track)
        self.assertFalse(args_track.choose)
        self.assertFalse(args_track.fixed)

    def test_parse_arguments_help_message(self):
        """Test help message generation."""
        with pytest.raises(SystemExit):
            parse_arguments(['--help'])
    
    def test_parse_arguments_custom_level(self):
        """Test custom pressure level setting."""
        test_args = ['input_file.nc', '--fixed', '--level', '500']
        
        args = parse_arguments(test_args)
        
        self.assertEqual(args.level, 500)
        self.assertTrue(args.fixed)
    
    def test_parse_arguments_track_options(self):
        """Test track vorticity and geopotential options."""
        test_args = [
            'input_file.nc', 
            '--track', 
            '--track_vorticity', 'max', 
            '--track_geopotential', 'max'
        ]
        
        args = parse_arguments(test_args)
        
        self.assertEqual(args.track_vorticity, 'max')
        self.assertEqual(args.track_geopotential, 'max')
        self.assertTrue(args.track)
    
    def test_parse_arguments_gfs_option(self):
        """Test GFS option."""
        test_args = ['input_file.nc', '--fixed', '--gfs']
        
        args = parse_arguments(test_args)
        
        self.assertTrue(args.gfs)
        self.assertTrue(args.fixed)
    
    def test_parse_arguments_outname_option(self):
        """Test custom output name."""
        test_args = ['input_file.nc', '--fixed', '--outname', 'custom_output']
        
        args = parse_arguments(test_args)
        
        self.assertEqual(args.outname, 'custom_output')
        self.assertTrue(args.fixed)
    
    def test_parse_arguments_save_nc_file_option(self):
        """Test save NetCDF file option."""
        test_args = ['input_file.nc', '--fixed', '--save_nc_file', 'False']
        
        args = parse_arguments(test_args)
        
        # Note: argparse with type=bool parses 'False' string as True because it's non-empty
        # This is a known argparse behavior
        self.assertTrue(args.save_nc_file)  # 'False' string evaluates to True
        
        # Test with default (no --save_nc_file argument)
        test_args_default = ['input_file.nc', '--fixed']
        args_default = parse_arguments(test_args_default)
        self.assertTrue(args_default.save_nc_file)  # default is True


class TestCliInterfaceIntegration(unittest.TestCase):
    """Integration tests for CLI interface functionality."""
    
    def test_module_imports(self):
        """Test that the module imports correctly."""
        import cli_interface
        self.assertTrue(hasattr(cli_interface, 'parse_arguments'))
    
    def test_argparse_integration(self):
        """Test argparse integration with various combinations."""
        # Test short form arguments
        args_short = parse_arguments(['test.nc', '-f'])
        self.assertTrue(args_short.fixed)
        
        args_short2 = parse_arguments(['test.nc', '-t'])  
        self.assertTrue(args_short2.track)
        
        args_short3 = parse_arguments(['test.nc', '-c'])
        self.assertTrue(args_short3.choose)
    
    def test_error_handling_patterns(self):
        """Test common error handling patterns."""
        # Test with missing required argument
        with pytest.raises(SystemExit):
            parse_arguments([])  # Missing infile
    
    def test_argument_types_and_validation(self):
        """Test argument types and validation."""
        # Test integer type for level
        args = parse_arguments(['input.nc', '--fixed', '--level', '925'])
        self.assertEqual(args.level, 925)
        self.assertIsInstance(args.level, int)
        
        # Test string type for filename
        args2 = parse_arguments(['my_file.nc', '--track'])
        self.assertEqual(args2.infile, 'my_file.nc')
        self.assertIsInstance(args2.infile, str)
        
        # Test boolean flags
        args3 = parse_arguments(['input.nc', '--verbose', '--fixed'])
        self.assertIsInstance(args3.verbose, bool)
        self.assertTrue(args3.verbose)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and special scenarios."""
    
    def test_special_filenames(self):
        """Test parsing with special filename characters."""
        special_files = [
            'file-with-dashes.nc',
            'file_with_underscores.nc',
            'file123numbers.nc',
            'UPPERCASE.NC'
        ]
        
        for filename in special_files:
            args = parse_arguments([filename, '--fixed'])
            self.assertEqual(args.infile, filename)
            self.assertTrue(args.fixed)
    
    def test_long_arguments_vs_short(self):
        """Test long vs short argument formats."""
        # Test long vs short forms
        args_long = parse_arguments(['input.nc', '--verbose', '--fixed'])
        args_short = parse_arguments(['input.nc', '-v', '-f'])
        
        self.assertEqual(args_long.verbose, args_short.verbose)
        self.assertEqual(args_long.fixed, args_short.fixed)
    
    def test_argument_order_independence(self):
        """Test that argument order doesn't matter for optional args."""
        # Test different orders of the same arguments
        args1 = parse_arguments(['input.nc', '--verbose', '--choose'])
        args2 = parse_arguments(['input.nc', '--choose', '--verbose'])
        
        self.assertEqual(args1.verbose, args2.verbose)
        self.assertEqual(args1.choose, args2.choose)
    
    def test_track_vorticity_choices(self):
        """Test track vorticity choice validation."""
        # Valid choices
        args1 = parse_arguments(['input.nc', '--track', '--track_vorticity', 'min'])
        self.assertEqual(args1.track_vorticity, 'min')
        
        args2 = parse_arguments(['input.nc', '--track', '--track_vorticity', 'max'])  
        self.assertEqual(args2.track_vorticity, 'max')
        
        # Invalid choice should cause SystemExit
        with pytest.raises(SystemExit):
            parse_arguments(['input.nc', '--track', '--track_vorticity', 'invalid'])
    
    def test_track_geopotential_choices(self):
        """Test track geopotential choice validation."""
        # Valid choices
        args1 = parse_arguments(['input.nc', '--track', '--track_geopotential', 'min'])
        self.assertEqual(args1.track_geopotential, 'min')
        
        args2 = parse_arguments(['input.nc', '--track', '--track_geopotential', 'max'])
        self.assertEqual(args2.track_geopotential, 'max')
        
        # Invalid choice should cause SystemExit
        with pytest.raises(SystemExit):
            parse_arguments(['input.nc', '--track', '--track_geopotential', 'invalid'])


if __name__ == '__main__':
    unittest.main()
