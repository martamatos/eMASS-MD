import os
import unittest

import pandas as pd

from src.parsing.add_ligand_to_mol2 import add_info_to_mol2_file
from src.parsing.change_atom_types_mol2 import change_atom_types_mol2
from src.parsing.parse_vmd_output import parse_vmd_rmsd_output
from src.parsing.reorder_ligand_atoms import reorder_ligand_atoms_from_off, rename_ligand_atoms
from src.parsing.substitute_ligands import substitute_ligands
from src.tests import get_tests_dir


def _get_file_content(filename):
    f_in = open(filename, 'r')
    res = f_in.read()
    f_in.close()

    return res


class TestAddInfoToMol2File(unittest.TestCase):
    
    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)
        

    def test_add_info_to_mol2_file(self):

        filename = 'test_files/test_parsing/add_ligand_to_mol2/TALB_WT_APO-001_S7P_flexible_scored'
        ligand = 'S7P'
        new_atom_list = ['H', 'C.3', 'C.3', 'C.3', 'O.3', 'P.3', 'O.co2', 'C.3', 'C.3', 'O.3',
                         'O.co2', 'O.co2', 'O.3', 'H', 'H', 'O.3', 'H', 'H', 'H', 'H', 'H', 'H',
                         'C.2', 'C.3', 'H', 'C.3', 'O.2', 'O.3', 'H', 'H', 'H']

        add_info_to_mol2_file(filename, ligand, new_atom_list)

        true_res = _get_file_content('test_files/test_parsing/add_ligand_to_mol2/TALB_WT_APO-001_S7P_flexible_scored2_true_res.mol2')
        res = _get_file_content('test_files/test_parsing/add_ligand_to_mol2/TALB_WT_APO-001_S7P_flexible_scored2.mol2')

        self.assertEqual(res, true_res)

    def test_add_info_to_mol2_file_no_atoms(self):

        filename = 'test_files/test_parsing/add_ligand_to_mol2/G6PD_WT_APO-001_G6P_flexible_scored'
        ligand = 'G6P'


        add_info_to_mol2_file(filename, ligand)

        true_res = _get_file_content('test_files/test_parsing/add_ligand_to_mol2/G6PD_WT_APO-001_G6P_flexible_scored2_true_res.mol2')
        res = _get_file_content('test_files/test_parsing/add_ligand_to_mol2/G6PD_WT_APO-001_G6P_flexible_scored2.mol2')

        self.assertEqual(res, true_res)


class TestChangeAtomTypesTest(unittest.TestCase):
    
    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)


    def test_change_atom_types_mol2(self):
        filename = '../tests/test_files/test_parsing/change_atom_types_mol2/TALB_WT_APO-400_S7P_flexible_scored'
        ligand = 'S7P'
        new_atom_list = ['H', 'C.3', 'C.3', 'C.3', 'O.3', 'P.3', 'O.co2', 'C.3', 'C.3', 'O.3',
                         'O.co2', 'O.co2', 'O.3', 'H', 'H', 'O.3', 'H', 'H', 'H', 'H', 'H', 'H',
                         'C.2', 'C.3', 'H', 'C.3', 'O.2', 'O.3', 'H', 'H', 'H']

        change_atom_types_mol2(filename, ligand, new_atom_list)

        true_res = _get_file_content('test_files/test_parsing/change_atom_types_mol2/TALB_WT_APO-400_S7P_flexible_scored2_true_res.mol2')
        res = _get_file_content('test_files/test_parsing/change_atom_types_mol2/TALB_WT_APO-400_S7P_flexible_scored2.mol2')

        self.assertEqual(res, true_res)


class TestParseVmdRmsdOutput(unittest.TestCase):
    
    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)

    def test_parse_vmd_rmsd_output(self):
        true_res = pd.read_csv('test_files/test_parsing/parse_vmd/true_res_trajrmsd_dcd_mod.dat', sep='\t')

        parse_vmd_rmsd_output('test_files/test_parsing/parse_vmd/trajrmsd_dcd.dat')
        res = pd.read_csv('test_files/test_parsing/parse_vmd/trajrmsd_dcd_mod.dat', sep='\t')

        self.assertEquals(true_res.values.tolist()[1:], res.values.tolist()[1:])


class TestReorderLigandAtomsTest(unittest.TestCase):
    
    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)
        

    def test_reorder_ligand_atoms_end(self):
        pdb_file = '../tests/test_files/test_parsing/reorder_ligand_atoms/TALB_WT_APO-033_receptor_F6P_docked'
        off_file = '../tests/test_files/test_parsing/reorder_ligand_atoms/F6P'
        ligand = 'F6P'
        n_atoms = 27

        reorder_ligand_atoms_from_off(pdb_file, off_file, ligand, n_atoms)

        true_res = _get_file_content('test_files/test_parsing/reorder_ligand_atoms/TALB_WT_APO-033_receptor_F6P_docked2_true_res.pdb')
        res = _get_file_content('test_files/test_parsing/reorder_ligand_atoms/TALB_WT_APO-033_receptor_F6P_docked2.pdb')

        self.assertEqual(res, true_res)

    def test_reorder_ligand_atoms_middle(self):
        pdb_file = '../tests/test_files/test_parsing/reorder_ligand_atoms/G6PD_WT_APO-001_receptor_6PGL_docked'
        off_file = '../tests/test_files/test_parsing/reorder_ligand_atoms/6PGL'
        ligand = '6PG'
        n_atoms = 25

        reorder_ligand_atoms_from_off(pdb_file, off_file, ligand, n_atoms)

        true_res = _get_file_content('test_files/test_parsing/reorder_ligand_atoms/G6PD_WT_APO-001_receptor_6PGL_docked2_true_res.pdb')
        res = _get_file_content('test_files/test_parsing/reorder_ligand_atoms/G6PD_WT_APO-001_receptor_6PGL_docked2.pdb')

        self.assertEqual(res, true_res)


    def test_reorder_ligand_atoms_middle_cofactor(self):
        pdb_file = '../tests/test_files/test_parsing/reorder_ligand_atoms/G6PD_WT_NADP-001_receptor_6PGL_docked'
        off_file = '../tests/test_files/test_parsing/reorder_ligand_atoms/6PGL'
        ligand = '6PG'
        n_atoms = 25

        reorder_ligand_atoms_from_off(pdb_file, off_file, ligand, n_atoms)

        true_res = _get_file_content('test_files/test_parsing/reorder_ligand_atoms/G6PD_WT_NADP-001_receptor_6PGL_docked2_true_res.pdb')
        res = _get_file_content('test_files/test_parsing/reorder_ligand_atoms/G6PD_WT_NADP-001_receptor_6PGL_docked2.pdb')

        self.assertEqual(res, true_res)


    def test_rename_ligand_atoms_auto(self):
        true_res = _get_file_content('test_files/test_parsing/rename_ligand_atoms/TALB_WT_APO-001_receptor_S7P_docked2_true_res.pdb')

        pdb_file = 'test_files/test_parsing/rename_ligand_atoms/TALB_WT_APO-001_receptor_S7P_docked'
        ligand_list = ['S7P']

        rename_ligand_atoms(pdb_file, ligand_list)

        res = _get_file_content('test_files/test_parsing/rename_ligand_atoms/TALB_WT_APO-001_receptor_S7P_docked2.pdb')
        self.assertEqual(res, true_res)


    def test_rename_ligand_atoms_from_map(self):
        true_res = _get_file_content('test_files/test_parsing/rename_ligand_atoms/G6PD_WT_APO-001_receptor_G6P_docked2_true_res.pdb')

        pdb_file = 'test_files/test_parsing/rename_ligand_atoms/G6PD_WT_APO-001_receptor_G6P_docked'
        ligand_list = ['G6P']
        map_atom_names = {'P1': 'P1', 'O7': 'O6', 'C16': 'C6', 'C13': 'C3', 'O2': 'O1', 'C15': 'C5', 'C14': 'C4',
                          'C12': 'C2', 'C11': 'C1', 'O6': 'O5', 'O5': 'O4', 'O4': 'O3', 'O3': 'O2', 'H21': 'H5',
                          'H20': 'H4', 'H17': 'H1', 'H19': 'H3', 'H18': 'H2', 'H27': 'H11', 'H26': 'H10', 'H25': 'H9',
                          'H24': 'H8', 'H22': 'H6', 'H23': 'H7', 'O8': 'O7', 'O9': 'O8', 'O10': 'O9'}

        rename_ligand_atoms(pdb_file, ligand_list, map_atom_names)

        res = _get_file_content('test_files/test_parsing/rename_ligand_atoms/G6PD_WT_APO-001_receptor_G6P_docked2.pdb')
        self.assertEqual(res, true_res)


class TestSubstituteLigandsTest(unittest.TestCase):

    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)

    def test_susbtitute_ligands(self):
        true_res = _get_file_content('test_files/test_parsing/susbtitute_ligands/GAPDH_WT_HALO_G3P_lasframe_noWATNa2_true_res.pdb')

        pdb_file = 'test_files/test_parsing/susbtitute_ligands/GAPDH_WT_HALO_G3P_lasframe_noWATNa'
        old_ligand = 'G3H'
        new_ligand = 'DPG'
        atom_map = {'P1': 'P2', 'O3': 'O8', 'O4': 'O9', 'O5': 'O10', 'C1': 'C3', 'C2': 'C2', 'C3': 'C1', 'H1': 'H3',
                    'H2': 'H4', 'O2': 'O6', 'H3': 'H1', 'H5': 'H2', 'O6': 'O5'}

        substitute_ligands(pdb_file, old_ligand, new_ligand, atom_map)

        res = _get_file_content('test_files/test_parsing/susbtitute_ligands/GAPDH_WT_HALO_G3P_lasframe_noWATNa2.pdb')
        self.assertEqual(res, true_res)
