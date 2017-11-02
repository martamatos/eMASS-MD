import os
import unittest

import pandas as pd

from src.kinetics_integration.binning import bin_keqs
from src.kinetics_integration.prepare_run_fit_script import prep_run_fit_script
from src.kinetics_integration.utils import import_rateconstants, set_df_column_labels
from src.tests import get_tests_dir


def _get_file_content(filename):
    with open(filename, 'r') as f_in:
        res = f_in.read()

    return res


"""
class TestPrepRunFitScript(unittest.TestCase):

    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)


    def test_prep_run_fit_script(self):
        true_res = _get_file_content('test_files/test_kinetics_integration/prep_run_fit_script/true_res_run_kinetic_fit.sh')

        main_dir = '/zhome/89/0/90554/kinetics_fit/'
        file_out = 'test_files/test_kinetics_integration/prep_run_fit_script/run_kinetic_fit.sh'
        enzyme = 'G6PDH2r'
        enzyme_folder = 'G6PD'
        job_number = 2
        dataset_list = ['Km', 'Keq', 'KmKeqdKd1', 'KmKeqdKd2']
        dKd_list = ['1.e-9', '1.e-8', '1.e-7', '0.000001', '0.00001', '0.0001', '0.001', '0.01', '0.1', '1.', '10', '100',
                    '1000', '10000', '100000', '1000000', '10000000', '100000000', '1000000000']
        n_trials = 20
        walltime = '5:00:00'

        prep_run_fit_script(file_out, main_dir, enzyme, enzyme_folder, python_packages_path, job_number, dataset_list,
                            dKd_list, n_trials, walltime=walltime)
            (file_out, main_dir, data_file, enzyme_folder, python_packages_path, job_number, n_trials,
                        start_fit, end_fit, walltime='5:00:00', qsub=False, fit_label='')
        res = _get_file_content(file_out)
        self.assertEqual(true_res.strip(), res.strip())
"""


class TestUtils(unittest.TestCase):

    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)

    def test_set_df_column_names(self):
        true_res = pd.read_csv('test_files/test_kinetics_integration/utils/true_res_set_df_column_names.csv', sep='\t', index_col=0)
        true_res = true_res.ix[:, 1:]

        data_df = pd.read_csv('test_files/test_kinetics_integration/utils/rateconst_ENO_dKd__all_0.25.dat', sep='\t', header=None)
        column_order = [0, 1, 4, 5, 2, 3]
        data_df = data_df.ix[:, 1:]
        data_df = set_df_column_labels(data_df, column_order=column_order)

        data_df.to_csv('test_files/test_kinetics_integration/utils/set_df_column_names_temp.csv', sep='\t')
        data_df = pd.read_csv('test_files/test_kinetics_integration/utils/set_df_column_names_temp.csv', sep='\t', index_col=0)

        self.assertTrue(true_res.equals(data_df))

    def test_import_data(self):
        true_res = pd.read_csv('test_files/test_kinetics_integration/utils/true_res_import_data.csv', sep='\t', index_col=0)
        file_in = 'test_files/test_kinetics_integration/utils/rateconst_ENO_dKd__all_0.25.dat'
        column_order = [0, 1, 4, 5, 2, 3]
        ssd_threshold = 1

        data_df, fitness = import_rateconstants(file_in, ssd_threshold, None, column_order=column_order)

        data_df.to_csv('test_files/test_kinetics_integration/utils/import_data_temp.csv', sep='\t', index_col=0)
        data_df = pd.read_csv('test_files/test_kinetics_integration/utils/import_data_temp.csv', sep='\t', index_col=0)

        self.assertTrue(true_res.equals(data_df))

    def test_import_data_as_ratios(self):
        true_res = pd.read_csv('test_files/test_kinetics_integration/utils/true_res_import_data_as_ratios.csv', sep='\t', index_col=0)
        file_in = 'test_files/test_kinetics_integration/utils/rateconst_ENO_dKd__all_0.25.dat'
        column_order = [0, 2, 1]
        ssd_threshold = 1

        data_df, fitness = import_rateconstants(file_in, ssd_threshold, None, column_order=column_order, convert_to_ratios=True)

        data_df.to_csv('test_files/test_kinetics_integration/utils/import_data_as_ratios_temp.csv', sep='\t', index_col=0)
        data_df = pd.read_csv('test_files/test_kinetics_integration/utils/import_data_as_ratios_temp.csv', sep='\t', index_col=0)

        self.assertTrue(true_res.equals(data_df))


class TestBinKeqs(unittest.TestCase):

    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)


    def test_bin_keqs_fixed_size(self):
        true_res_pattern_df = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/true_res_bin_keqs_pattern_df_fixed_size.csv', sep='\t')
        true_res_pattern_df_by_model = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/true_res_bin_keqs_pattern_df_by_model_fixed_size.csv', sep='\t')

        column_order = [0, 2, 1]
        convert_to_ratios = True
        ssd_threshold = 1
        scale_data = False
        n_samples_per_bin = 5

        file_in = 'test_files/test_kinetics_integration/bin_keqs/rateconst_ENO_dKd_1.csv'
        file_out = 'test_files/test_kinetics_integration/bin_keqs/ENO_'

        data_df, fitness = import_rateconstants(file_in, ssd_threshold, None, column_order,
                                                convert_to_ratios=convert_to_ratios)
        pattern_df, pattern_df_by_model = bin_keqs(data_df, fixed_size=True, n_samples_per_bin=n_samples_per_bin,
                                                   file_out=file_out)

        pattern_df.to_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_fixed_size.csv', sep='\t')
        pattern_df_by_model.to_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_by_model_fixed_size.csv', sep='\t')

        pattern_df = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_fixed_size.csv', sep='\t')
        pattern_df_by_model = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_by_model_fixed_size.csv', sep='\t')

        self.assertTrue(true_res_pattern_df.equals(pattern_df))
        self.assertTrue(true_res_pattern_df_by_model.equals(pattern_df_by_model))


    def test_bin_keqs_fixed_size_tala2(self):
        true_res_pattern_df = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/true_res_bin_keqs_pattern_df_fixed_size_tala2.csv', sep='\t')
        true_res_pattern_df_by_model = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/true_res_bin_keqs_pattern_df_by_model_fixed_size_tala2.csv', sep='\t')

        column_order = [0, 3, 5, 2, 4, 1]
        convert_to_ratios = True
        ssd_threshold = 1
        scale_data = False
        n_samples_per_bin = 5

        file_in = 'test_files/test_kinetics_integration/bin_keqs/rateconst_TALA2_nodKd_1.csv'
        file_out = 'test_files/test_kinetics_integration/bin_keqs/TALA2_'

        data_df, fitness = import_rateconstants(file_in, ssd_threshold, None, column_order,
                                                convert_to_ratios=convert_to_ratios)
        pattern_df, pattern_df_by_model = bin_keqs(data_df, fixed_size=True, n_samples_per_bin=n_samples_per_bin,
                                                   file_out=file_out)

        pattern_df.to_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_fixed_size_tala2.csv', sep='\t')
        pattern_df_by_model.to_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_by_model_fixed_size_tala2.csv', sep='\t')

        pattern_df = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_fixed_size_tala2.csv', sep='\t')
        pattern_df_by_model = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_by_model_fixed_size_tala2.csv', sep='\t')

        self.assertTrue(true_res_pattern_df.equals(pattern_df))
        self.assertTrue(true_res_pattern_df_by_model.equals(pattern_df_by_model))


    def test_bin_keqs_fixed_size_tala22(self):
        true_res_pattern_df = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/true_res_bin_keqs_pattern_df_fixed_size_tala22.csv', sep='\t')
        true_res_pattern_df_by_model = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/true_res_bin_keqs_pattern_df_by_model_fixed_size_tala22.csv', sep='\t')

        column_order = [0, 3, 5, 2, 4, 1]
        convert_to_ratios = True
        ssd_threshold = 1
        scale_data = False
        n_samples_per_bin = 5

        file_in = 'test_files/test_kinetics_integration/bin_keqs/rateconst_TALA2_dKd_1.csv'
        file_out = 'test_files/test_kinetics_integration/bin_keqs/TALA2_'

        data_df, fitness = import_rateconstants(file_in, ssd_threshold, None, column_order,
                                                convert_to_ratios=convert_to_ratios)

        pattern_df, pattern_df_by_model = bin_keqs(data_df, fixed_size=True, n_samples_per_bin=n_samples_per_bin,
                                                   file_out=file_out)


        pattern_df.to_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_fixed_size_tala22.csv', sep='\t')
        pattern_df_by_model.to_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_by_model_fixed_size_tala22.csv', sep='\t')

        pattern_df = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_fixed_size_tala22.csv', sep='\t')
        pattern_df_by_model = pd.read_csv('test_files/test_kinetics_integration/bin_keqs/bin_keqs_pattern_df_by_model_fixed_size_tala22.csv', sep='\t')

        self.assertTrue(true_res_pattern_df.equals(pattern_df))
        self.assertTrue(true_res_pattern_df_by_model.equals(pattern_df_by_model))