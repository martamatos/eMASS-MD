import os
import unittest

import pandas as pd
import shutil
from src.process_traj.analyze_traj import dist_backbone_ligand, get_atoms_dist, write_vmd_script, process_out_files, \
    process_out_files_min
from src.process_traj.rmsd_MDAnalysis import calculate_traj_rmsd_fit, average_rmsd_from_file
from src.process_traj.traj_utils import convert_rst_to_pdb, weed_traj
from src.tests import get_tests_dir


def _get_file_content(filename):
    f_in = open(filename, "r")
    res = f_in.read()
    f_in.close()

    return res
    
    
class TestAnalyzeTraj(unittest.TestCase):

    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)

    def process_out_files(self):
        return 0

    def test_dist_backbone_ligand_compare_vmd(self):
        true_res = []
        for i in range(6):
            df = pd.read_csv(''.join(['test_files/test_process_traj/analyze_traj/true_bond', str(i+1), '.dat']), sep='\t', header=None)
            true_res.append(df[1].values.tolist())

        traj = 'test_files/G6PD_WT_APO_G6P_docked_cl111_10frames.dcd'
        prmtop = 'test_files/G6PD_WT_APO_G6P_docked.prmtop'
        ligand = 'G6P'
        file_out = 'test_files/test_process_traj/analyze_traj/bond'

        dist_backbone_ligand(traj, prmtop, ligand, file_out)

        df = pd.read_csv('.'.join([file_out, 'dat']), sep='\t')
        del df['frame']

        for i in range(6):
            for j in range(len(true_res[i])):
                self.assertAlmostEquals(true_res[i][j], df[df.columns[i]].values.tolist()[j], 5)

    def test_dist_backbone_ligand_relative(self):
        true_res = pd.read_csv(''.join(['test_files/test_process_traj/analyze_traj/true_res_bond_relative.dat']), sep='\t')

        traj = 'test_files/G6PD_WT_APO_G6P_docked_cl111_10frames.dcd'
        prmtop = 'test_files/G6PD_WT_APO_G6P_docked.prmtop'
        ligand = 'G6P'
        file_out = 'test_files/test_process_traj/analyze_traj/bond'

        dist_backbone_ligand(traj, prmtop, ligand, file_out)

        df = pd.read_csv('_'.join([file_out, 'relative.dat']), sep='\t')

        self.assertListEqual(true_res.columns.values.tolist(), df.columns.values.tolist())
        for i in range(len(true_res.columns.values)):
            self.assertEquals(true_res[true_res.columns[i]].values.tolist(), df[df.columns[i]].values.tolist())

    def test_get_atom_dist_vmd(self):
        true_res = []
        for i in range(6):
            df = pd.read_csv(''.join(['test_files/test_process_traj/analyze_traj/true_bond', str(i+1), '.dat']), sep='\t', header=None)
            true_res.append(df[1].values.tolist())

        traj = 'test_files/G6PD_WT_APO_G6P_docked_cl111_10frames.dcd'
        prmtop = 'test_files/G6PD_WT_APO_G6P_docked.prmtop'

        # get atoms distance
        atom_pairs = [(7649, 2362), (7650, 2405), (7652, 2416), (7662, 6998)]
        file_out = 'test_files/test_process_traj/analyze_traj/MG_restraints'
        get_atoms_dist(traj, prmtop, atom_pairs, file_out)

        df = pd.read_csv('.'.join([file_out, 'dat']), sep='\t')
        del df['frame']

        for i in range(4):
            for j in range(len(true_res[i])):
                self.assertAlmostEquals(true_res[i][j], df[df.columns[i]].values.tolist()[j], 5)

    def test_get_atom_dist_relative(self):
        true_res = pd.read_csv(''.join(['test_files/test_process_traj/analyze_traj/true_res_MG_restraints_relative.dat']), sep='\t')

        traj = 'test_files/G6PD_WT_APO_G6P_docked_cl111_10frames.dcd'
        prmtop = 'test_files/G6PD_WT_APO_G6P_docked.prmtop'

        atom_pairs = [(7649, 2362), (7650, 2405), (7652, 2416), (7662, 6998)]
        file_out = 'test_files/test_process_traj/analyze_traj/MG_restraints'
        get_atoms_dist(traj, prmtop, atom_pairs, file_out)

        df = pd.read_csv('_'.join([file_out, 'relative.dat']), sep='\t')

        self.assertListEqual(['frame', '(7649, 2362)', '(7650, 2405)', '(7652, 2416)', '(7662, 6998)'], df.columns.values.tolist())
        for i in range(len(true_res.columns.values)):
            self.assertEquals(true_res[true_res.columns[i]].values.tolist(), df[df.columns[i]].values.tolist())


    def test_process_out_files_min(self):
        md_dir = 'test_files/test_process_traj/analyze_traj/process_out_files/'

        data_df_true_res = pd.read_csv(''.join([self.base_dir, '/', md_dir, 'true_res_data_df_min.csv']), sep='\t')
        fluctuations_df_true_res = pd.read_csv(''.join([self.base_dir, '/', md_dir, 'true_res_fluctuations_df_min.csv']), sep='\t')

        file_in = '/'.join([self.base_dir, md_dir, 'G6PD_WT_NADP_G6P_min_v0_l1.out'])
        file_out = 'data/min'
        data_dic, fluctuations_dic = process_out_files_min(md_dir, file_in, file_out)

        data_df = pd.DataFrame.from_dict(data_dic)
        fluctuations_df = pd.DataFrame.from_dict(fluctuations_dic, orient='index')

        data_df.to_csv(''.join([self.base_dir, '/', md_dir, '/data/data_df_min.csv']), sep='\t')
        fluctuations_df.to_csv(''.join([self.base_dir, '/', md_dir, '/data/fluctuations_df_min.csv']), sep='\t')

        data_df = pd.read_csv(''.join([self.base_dir, '/', md_dir, '/data/data_df_min.csv']), sep='\t')
        fluctuations_df = pd.read_csv(''.join([self.base_dir, '/', md_dir, '/data/fluctuations_df_min.csv']), sep='\t')

        self.assertTrue(data_df_true_res.equals(data_df))
        self.assertTrue(fluctuations_df_true_res.equals(fluctuations_df))

        if os.path.isdir('data'):
            shutil.rmtree('data')


    def test_process_out_files_equil(self):
        md_dir = 'test_files/test_process_traj/analyze_traj/process_out_files/'

        data_df_true_res = pd.read_csv(''.join([self.base_dir, '/', md_dir, 'true_res_data_df_equil.csv']), sep='\t')
        fluctuations_df_true_res = pd.read_csv(''.join([self.base_dir, '/', md_dir, 'true_res_fluctuations_df_equil.csv']), sep='\t')

        file_in = '/'.join([self.base_dir, md_dir, 'G6PD_WT_NADP_G6P_equil_v0_l1.out'])
        file_out = 'data/equil'
        data_dic, fluctuations_dic = process_out_files(md_dir, file_in, file_out)

        data_df = pd.DataFrame.from_dict(data_dic)
        fluctuations_df = pd.DataFrame.from_dict(fluctuations_dic, orient='index')

        data_df.to_csv(''.join([self.base_dir, '/', md_dir, '/data/data_df_equil.csv']), sep='\t')
        fluctuations_df.to_csv(''.join([self.base_dir, '/', md_dir, '/data/fluctuations_df_equil.csv']), sep='\t')

        data_df = pd.read_csv(''.join([self.base_dir, '/', md_dir, '/data/data_df_equil.csv']), sep='\t')
        fluctuations_df = pd.read_csv(''.join([self.base_dir, '/', md_dir, '/data/fluctuations_df_equil.csv']), sep='\t')

        self.assertTrue(data_df_true_res.equals(data_df))
        self.assertTrue(fluctuations_df_true_res.equals(fluctuations_df))

        if os.path.isdir('data'):
            shutil.rmtree('data')


    def test_process_out_files_prod(self):
        md_dir = 'test_files/test_process_traj/analyze_traj/process_out_files/'

        data_df_true_res = pd.read_csv(''.join([self.base_dir, '/', md_dir, 'true_res_data_df_prod.csv']), sep='\t')
        fluctuations_df_true_res = pd.read_csv(''.join([self.base_dir, '/', md_dir, 'true_res_fluctuations_df_prod.csv']), sep='\t')

        file_in = '/'.join([self.base_dir, md_dir, 'G6PD_WT_NADP_G6P_prod_v0_l1.out'])
        file_out = 'data/prod'
        data_dic, fluctuations_dic = process_out_files(md_dir, file_in, file_out)

        data_df = pd.DataFrame.from_dict(data_dic)
        fluctuations_df = pd.DataFrame.from_dict(fluctuations_dic, orient='index')

        data_df.to_csv(''.join([self.base_dir, '/', md_dir, '/data/data_df_prod.csv']), sep='\t')
        fluctuations_df.to_csv(''.join([self.base_dir, '/', md_dir, '/data/fluctuations_df_prod.csv']), sep='\t')

        data_df = pd.read_csv(''.join([self.base_dir, '/', md_dir, '/data/data_df_prod.csv']), sep='\t')
        fluctuations_df = pd.read_csv(''.join([self.base_dir, '/', md_dir, '/data/fluctuations_df_prod.csv']), sep='\t')

        self.assertTrue(data_df_true_res.equals(data_df))
        self.assertTrue(fluctuations_df_true_res.equals(fluctuations_df))

        if os.path.isdir('data'):
            shutil.rmtree('data')




    def test_write_vmd_script(self):
        true_res = _get_file_content('test_files/test_process_traj/analyze_traj/true_res_load_traj_no_label.tcl')

        traj = 'test_files/G6PD_WT_APO_G6P_docked_cl111_10frames.dcd'
        prmtop = 'test_files/G6PD_WT_APO_G6P_docked.prmtop'
        file_out = 'test_files/test_process_traj/analyze_traj/load_traj_no_label.tcl'
        ligand = 'G6P'
        write_vmd_script(traj, prmtop, file_out, ligand, atom_pairs=None)
        res = _get_file_content('test_files/test_process_traj/analyze_traj/load_traj_no_label.tcl')

        self.assertEquals(true_res.strip(), res.strip())


    def test_write_vmd_script_labels(self):
        true_res = _get_file_content('test_files/test_process_traj/analyze_traj/true_res_load_traj.tcl')

        traj = 'test_files/G6PD_WT_APO_G6P_docked_cl111_10frames.dcd'
        prmtop = 'test_files/G6PD_WT_APO_G6P_docked.prmtop'
        file_out = 'test_files/test_process_traj/analyze_traj/load_traj.tcl'
        atom_pairs = [(7648, 2361), (7649, 2404), (7651, 2415), (7661, 6997), (7665, 7048), (7669, 7064)]
        ligand = 'G6P'
        write_vmd_script(traj, prmtop, file_out, ligand, atom_pairs)
        res = _get_file_content('test_files/test_process_traj/analyze_traj/load_traj.tcl')

        self.assertEquals(true_res.strip(), res.strip())


class TestCalculateRMSDMDAnalysis(unittest.TestCase):

    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)


    def test_calculate_traj_rmsd_backbone(self):
        true_res = pd.read_csv('test_files/test_process_traj/calculate_rmsd_MDanalysis/true_res_trajrmsd_dcd.dat', sep='\t')
        true_res = true_res['mol0'].values.tolist()

        traj = 'test_files/G6PD_WT_APO_G6P_docked_cl111_10frames.dcd'
        prmtop = 'test_files/G6PD_WT_APO_G6P_docked.prmtop'
        file_out = 'test_files/test_process_traj/calculate_rmsd_MDanalysis/rmsdfit'
        selection = 'backbone'
        calculate_traj_rmsd_fit(traj, prmtop, file_out, selection)

        res = pd.read_csv('test_files/test_process_traj/calculate_rmsd_MDanalysis/rmsdfit.dat', sep=' ', header=None)
        res = res[res.columns[2]].values.tolist()

        for i in range(1, len(res)):
            self.assertAlmostEquals(round(float(res[i]), 3), round(float(true_res[i]), 3), delta=0.07)


    def test_calculate_traj_rmsd_ligand(self):
        true_res = pd.read_csv('test_files/test_process_traj/calculate_rmsd_MDanalysis/true_res_trajrmsd_ligand.dat', sep='\t')
        true_res = true_res['mol0'].values.tolist()

        traj = 'test_files/G6PD_WT_APO_G6P_docked_cl111_10frames.dcd'
        prmtop = 'test_files/G6PD_WT_APO_G6P_docked.prmtop'
        file_out = 'test_files/test_process_traj/calculate_rmsd_MDanalysis/rmsdfit_ligand'
        selection = 'resname G6P'
        calculate_traj_rmsd_fit(traj, prmtop, file_out, selection)

        res = pd.read_csv('test_files/test_process_traj/calculate_rmsd_MDanalysis/rmsdfit_ligand.dat', sep=' ', header=None)
        res = res[res.columns[2]].values.tolist()

        for i in range(1, len(res)):
            self.assertAlmostEquals(round(float(res[i]), 3), round(float(true_res[i]), 3), delta=0.3)


    def test_compare_vmd_res(self):
        vmd_res_crd = _get_file_content('test_files/test_process_traj/calculate_rmsd_MDanalysis/trajrmsd_crd.dat')
        vmd_res_dcd = _get_file_content('test_files/test_process_traj/calculate_rmsd_MDanalysis/trajrmsd_dcd.dat')

        self.assertEquals(vmd_res_crd.strip(), vmd_res_dcd.strip())


    def test_average_rmsd_from_file(self):
        file_rmsd = 'test_files/test_process_traj/calculate_rmsd_MDanalysis/rmsdfit.dat'
        file_out = 'test_files/test_process_traj/calculate_rmsd_MDanalysis/average_rmsd'
        average_rmsd_from_file(file_rmsd, file_out)


class TestTrajUtils(unittest.TestCase):

    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)

    # IMPORTANT!!!
    #  to run this test you need to comment all lines after "# go to the folder where the md simulation is"
    #  and the line 'import pytraj as pt'
    def test_weed_traj(self):
        true_res = _get_file_content('test_files/test_process_traj/process_traj/true_res_weed_traj_100')

        md_traj_dir = 'test_files'
        script_dir = 'test_files/test_process_traj/process_traj'

        base_name = 'G6PD_WT_APO_G6P_docked'
        md_sim = 1
        cl = 'cl100'
        frame_intervals = [100]

        weed_traj(md_traj_dir, script_dir, base_name, frame_intervals, md_sim, cl)
        res = _get_file_content('/'.join([self.base_dir, 'test_files/test_process_traj/process_traj/weed_traj.in']))

        self.assertEqual(true_res.strip(), res.strip())

    def test_convert_rst_to_pdb(self):
        true_res = 'test_files/test_process_traj/traj_utils/convert_rst_to_pdb/true_res_cpp_traj.in'

        rst_file = 'test_files/test_process_traj/traj_utils/convert_rst_to_pdb/G6PD_WT_NADP_G6P_min_v0_l1.rst'
        prmtop = 'test_files/test_process_traj/traj_utils/convert_rst_to_pdb/G6PD_WT_NADP_G6P.prmtop'
        cpp_traj_file = 'test_files/test_process_traj/traj_utils/convert_rst_to_pdb/cpp_traj.in'
        convert_rst_to_pdb(rst_file, prmtop, cpp_traj_file=cpp_traj_file)

        with open(true_res, 'r') as f_in:
            true_res_content = f_in.read()

        with open(cpp_traj_file, 'r') as f_in:
            res_content = f_in.read()

        self.assertEquals(true_res_content, res_content)

