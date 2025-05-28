import pytest
from unittest.mock import MagicMock, patch
import pandas as pd
import os

from src.output_management import save_results_csv

def test_save_results_csv_success():
    results_subdirectory = '/fake/dir'
    app_logger = MagicMock()

    # Criando DataFrames mockados (pode ser dataframe real vazio)
    df1 = pd.DataFrame({'a': [1,2]})
    df2 = {'key1': [1,2]}  # Para testar o termo 'ResQ' que converte dict para df
    df3 = pd.DataFrame({'b': [3,4]})

    results_df_dictionary = {
        'AdvHTemp': df1,
        'ResQ': df2,
        'Zeta': df3,
        'UnknownTerm': df1
    }

    # Mock de os.path.join
    def mock_join(*args):
        return '/'.join(args)

    with patch('pandas.DataFrame.to_csv') as mock_to_csv, \
        patch('os.path.join', mock_join):  # Mock corrigido de os.path.join

        save_results_csv(results_df_dictionary, results_subdirectory, app_logger)

        # Agora imprime as chamadas do logger, para depurar:
        print("Logger info calls:")
        for call in app_logger.info.call_args_list:
            print(call.args[0])

        # Verificar se to_csv foi chamado 4 vezes (uma para cada termo)
        assert mock_to_csv.call_count == 4

        # Verificar se o logger info foi chamado com os caminhos corretos
        expected_calls = [
            (f'{results_subdirectory}/heat_terms/AdvHTemp.csv',),
            (f'{results_subdirectory}/ResQ.csv',),
            (f'{results_subdirectory}/vorticity_terms/Zeta.csv',),
            (f'{results_subdirectory}/UnknownTerm.csv',),
        ]

        info_calls = [call.args for call in app_logger.info.call_args_list]

        # Verificar que o caminho esperado está sendo incluído nas chamadas info
        for expected in expected_calls:
            print(f"Checking: {expected[0]}")  # Para depuração
            assert any(expected[0] in message for (message,) in info_calls)

        # Verificar se o logger error foi chamado para 'UnknownTerm'
        app_logger.error.assert_any_call("Unknown term: UnknownTerm. Saving to main directory.")

def test_save_results_csv_raises_and_logs_error():
    results_subdirectory = '/fake/dir'
    app_logger = MagicMock()

    df = pd.DataFrame({'a': [1,2]})
    results_df_dictionary = {'AdvHTemp': df}

    with patch('pandas.DataFrame.to_csv', side_effect=Exception("Disk full")), \
         patch('os.path.join', side_effect=lambda *args: '/'.join(args)):

        save_results_csv(results_df_dictionary, results_subdirectory, app_logger)

        # Verificar se logger error foi chamado com mensagem de erro
        assert any("Error saving CSV files:" in call.args[0] for call in app_logger.error.call_args_list)
