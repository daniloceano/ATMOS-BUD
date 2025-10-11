import os
import pytest
from unittest.mock import patch, MagicMock
from collections import namedtuple

from src.output_management import manage_output

# Definir uma classe simples para simular args
Args = namedtuple('Args', ['outname', 'fixed', 'track'])

@patch('src.output_management.backup_file')
@patch('os.makedirs')
def test_manage_output_basic(mock_makedirs, mock_backup):
    # Caso 1: args.outname não fornecido, fixed=True
    args = Args(outname=None, fixed=True, track=False)
    infile = 'data/input_file.nc'
    method = 'fixed'

    results_subdir, figs_subdir, outfile_name = manage_output(args, infile, method)

    # Nome esperado do arquivo
    expected_outfile = 'input_file_fixed'

    # Testar os retornos
    assert results_subdir == os.path.join('./Results/', expected_outfile)
    assert figs_subdir == os.path.join(results_subdir, 'Figures')
    assert outfile_name == expected_outfile

    # Verificar se os diretórios foram tentados criar com os caminhos corretos
    expected_dirs = [
        './Results/',
        results_subdir,
        figs_subdir,
        os.path.join(results_subdir, 'heat_terms'),
        os.path.join(results_subdir, 'moisture_terms'),
        os.path.join(results_subdir, 'vorticity_terms'),
    ]

    calls = [((d,), {'exist_ok': True}) for d in expected_dirs]
    mock_makedirs.assert_has_calls(calls, any_order=True)

    # Verificar se backup_file foi chamado com 'inputs/box_limits' e results_subdir
    mock_backup.assert_called_once_with('inputs/box_limits', results_subdir)


@patch('src.output_management.backup_file')
@patch('os.makedirs')
def test_manage_output_with_outname_and_track(mock_makedirs, mock_backup):
    # Caso 2: args.outname fornecido, track=True
    args = Args(outname='custom_output', fixed=False, track=True)
    infile = 'data/ignored.nc'
    method = 'methodB'

    results_subdir, figs_subdir, outfile_name = manage_output(args, infile, method)

    assert results_subdir == os.path.join('./Results/', 'custom_output')
    assert figs_subdir == os.path.join(results_subdir, 'Figures')
    assert outfile_name == 'custom_output'

    # Verificar se backup_file foi chamado com 'inputs/track' e results_subdir
    mock_backup.assert_called_once_with('inputs/track', results_subdir)

