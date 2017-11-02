import os
import unittest

from src.cluster_dock.analyze_clusters import get_rep_frame
from src.cluster_dock.cluster_dock import cluster_docking_poses
from src.tests import get_tests_dir


def _get_file_content( filename):
    f_in = open(filename, "r")
    res = f_in.read()
    f_in.close()

    return res


class TestClusterDock(unittest.TestCase):

    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)

    def test_dock_cluster_manual_with_ref(self):
        true_top_clusters = ['cl102', 'cl112', 'cl111']
        true_distances = _get_file_content('test_files/test_cluster_dock/true_distances_wref')
        true_log = _get_file_content('test_files/test_cluster_dock/true_log_clustering_wref.txt')
        true_vmd_viz = _get_file_content('test_files/test_cluster_dock/true_viz_clusters_wref.tcl')
        true_cl1_list = _get_file_content('test_files/test_cluster_dock/true_cl102_frames_wref.txt')
        true_cl2_list = _get_file_content('test_files/test_cluster_dock/true_cl112_frames_wref.txt')
        true_cl3_list = _get_file_content('test_files/test_cluster_dock/true_cl111_frames_wref.txt')

        dock_base_folder = 'test_files/test_cluster_dock'

        chain = 'A'
        chain_ref = 'A'
        protein = 'GAPDH_WT_APO'
        ligand = 'G3P'
        ligand_ref = 'G3H'
        ligand_resname = 'H_G3H'

        ligand_resnum = 331
        ligand_resnum_ref = 350
        n_pdbs_docked = 50  # consider just counting the number of files in the directory
        ref_pdb = 'test_files/test_cluster_dock/1DC4_chainA.pdb'

        frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

        atom_pairs = [('O1', (209, 'CA')), ('C1', (176, 'CE1')), ('C2', (149, 'SG'))]
        atom_pairs_ref = [('O1P', (209, 'CA')), ('C3', (176, 'CE1')), ('C2', (149, 'SG'))]
        map_ligand_atoms = {'O1P': 'O1', 'C3': 'C1', 'C2': 'C2'}
        bandwidth_quantile = 0.2
        bond_cutoff = 12

        top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                              bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                              ligand_resname=ligand_resname, ligand_ref=ligand_ref,ligand_resname_ref=None,
                                              ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                              map_ligand_atoms=map_ligand_atoms)

        distances = _get_file_content('test_files/test_cluster_dock/distances')
        log = _get_file_content('test_files/test_cluster_dock/log_clustering.txt')
        vmd_viz = _get_file_content('test_files/test_cluster_dock/viz_clusters.tcl')
        cl1_list = _get_file_content('test_files/test_cluster_dock/cl102_frames.txt')
        cl2_list = _get_file_content('test_files/test_cluster_dock/cl112_frames.txt')
        cl3_list = _get_file_content('test_files/test_cluster_dock/cl111_frames.txt')

        self.assertListEqual(true_top_clusters, list(top_clusters))
        self.assertEquals(true_distances.strip(), distances.strip())
        self.assertEquals(true_log.strip(), log.strip())
        #self.assertEquals(true_vmd_viz.strip(), vmd_viz.strip())
        self.assertEquals(true_cl1_list.strip(), cl1_list.strip())
        self.assertEquals(true_cl2_list.strip(), cl2_list.strip())
        self.assertEquals(true_cl3_list.strip(), cl3_list.strip())


    def test_dock_cluster_manual(self):
        true_top_clusters = ['cl101', 'cl111', 'cl100']
        true_distances = _get_file_content('test_files/test_cluster_dock/true_distances')
        true_log = _get_file_content('test_files/test_cluster_dock/true_log_clustering.txt')
        true_vmd_viz = _get_file_content('test_files/test_cluster_dock/true_viz_clusters.tcl')
        true_cl1_list = _get_file_content('test_files/test_cluster_dock/true_cl101_frames.txt')
        true_cl2_list = _get_file_content('test_files/test_cluster_dock/true_cl111_frames.txt')
        true_cl3_list = _get_file_content('test_files/test_cluster_dock/true_cl100_frames.txt')

        dock_base_folder = 'test_files/test_cluster_dock'

        chain = 'A'
        chain_ref = 'A'
        protein = 'GAPDH_WT_APO'
        ligand = 'G3P'
        ligand_ref = 'G3H'
        ligand_resname = 'H_G3H'

        ligand_resnum = 331
        ligand_resnum_ref = 331
        n_pdbs_docked = 50  # consider just counting the number of files in the directory
        ref_pdb = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

        frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

        atom_pairs = [('O1', (209, 'CA')), ('C1', (176, 'CE1')), ('C2', (149, 'SG'))]
        atom_pairs_ref = [('O1', (209, 'CA')), ('C1', (176, 'CE1')), ('C2', (149, 'SG'))]
        map_ligand_atoms = {'O1': 'O1', 'C1': 'C1', 'C2': 'C2'}
        bandwidth_quantile = 0.2
        bond_cutoff = 12

        top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                              bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                              ligand_resname=ligand_resname, ligand_ref=ligand_ref,ligand_resname_ref=None,
                                              ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                              map_ligand_atoms=map_ligand_atoms)

        distances = _get_file_content('test_files/test_cluster_dock/distances')
        log = _get_file_content('test_files/test_cluster_dock/log_clustering.txt')
        vmd_viz = _get_file_content('test_files/test_cluster_dock/viz_clusters.tcl')
        cl1_list = _get_file_content('test_files/test_cluster_dock/cl101_frames.txt')
        cl2_list = _get_file_content('test_files/test_cluster_dock/cl111_frames.txt')
        cl3_list = _get_file_content('test_files/test_cluster_dock/cl100_frames.txt')

        self.assertListEqual(true_top_clusters, list(top_clusters))
        self.assertEquals(true_distances.strip(), distances.strip())
        self.assertEquals(true_log.strip(), log.strip())
        self.assertEquals(true_vmd_viz.strip(), vmd_viz.strip())
        self.assertEquals(true_cl1_list.strip(), cl1_list.strip())
        self.assertEquals(true_cl2_list.strip(), cl2_list.strip())
        self.assertEquals(true_cl3_list.strip(), cl3_list.strip())


class TestAnalyzeClusters(unittest.TestCase):

    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)

    def test_get_rep_frame(self):
        true_rep_frame = _get_file_content('test_files/test_cluster_dock/true_rep_frame.dat')
        
        dock_base_folder = 'test_files/test_cluster_dock'
        clusters = ['cl100_get_rep', 'cl101_get_rep', 'cl111_get_rep']
        atom_pairs = [('O1', (209, 'CA')), ('C1', (176, 'CE1')), ('C2', (149, 'SG'))]
        ligand = 'G3H'
        file_out = 'rep_frame'
        get_rep_frame(dock_base_folder, clusters, atom_pairs, ligand, file_out)

        rep_frame = _get_file_content('test_files/test_cluster_dock/rep_frame.dat')
        self.assertEquals(true_rep_frame.strip(), rep_frame.strip())
