import os
import unittest

from src.mmpbsa.harvest_mmpbsa_results import harvest_mmpbsa_perframe, average_dG, get_ddG
from src.mmpbsa.prepare_mmpbsa import prep_cpptraj_files, prep_get_frame_from_traj, prep_complex_pdbs, prep_mmpbsa_run, \
    prep_mmpbsa_in, pb_perframe
from src.tests import get_tests_dir


def _get_file_content(filename):
    with open(filename, 'r') as f_in:
        res = f_in.read()

    return res


class TestPrepareMMPBSA(unittest.TestCase):

    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)


    def test_prep_cpptraj_files_complex_solv(self):
        true_res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/true_res_get_last_traj_frame_complex_solv.in')

        output_dir = 'test_files/test_mmpbsa/prepare_mmpbsa'
        traj_name = 'G6PD_WT_NADP_6PGL_docked_cl111_100.crd'
        ligand = '6PG'
        cpp_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']

        prep_cpptraj_files(output_dir, traj_name, ligand, cpp_names)
        res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/get_traj_frame_complex_solv.in')

        self.assertEqual(true_res.strip(), res.strip())

    def test_prep_cpptraj_files_complex_vac(self):
        true_res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/true_res_get_last_traj_frame_complex_vac.in')

        output_dir = 'test_files/test_mmpbsa/prepare_mmpbsa/'
        traj_name = 'G6PD_WT_NADP_6PGL_docked_cl111_100.crd'
        ligand = '6PG'
        cpp_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']

        prep_cpptraj_files(output_dir, traj_name, ligand, cpp_names)
        res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/get_traj_frame_complex_vac.in')

        self.assertEqual(true_res.strip(), res.strip())

    def test_prep_cpptraj_files_ligand_vac(self):
        true_res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/true_res_get_last_traj_frame_ligand_vac.in')

        output_dir = 'test_files/test_mmpbsa/prepare_mmpbsa/'
        traj_name = 'G6PD_WT_NADP_6PGL_docked_cl111_100.crd'
        ligand = '6PG'
        cpp_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']

        prep_cpptraj_files(output_dir, traj_name, ligand, cpp_names)
        res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/get_traj_frame_ligand_vac.in')

        self.assertEqual(true_res.strip(), res.strip())

    def test_prep_cpptraj_files_receptor_vac(self):
        true_res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/true_res_get_last_traj_frame_receptor_vac.in')

        output_dir = 'test_files/test_mmpbsa/prepare_mmpbsa/'
        traj_name = 'G6PD_WT_NADP_6PGL_docked_cl111_100.crd'
        ligand = '6PG'
        cpp_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']

        prep_cpptraj_files(output_dir, traj_name, ligand, cpp_names)
        res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/get_traj_frame_receptor_vac.in')

        self.assertEqual(true_res.strip(), res.strip())

    def test_prep_get_frame_from_traj(self):
        true_res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/true_res_run_cpptraj_mmpbsa.sh')

        output_dir = 'test_files/test_mmpbsa/prepare_mmpbsa'
        md_traj_dir = '/home/marta/G6PD/3_MD_post_dock/2_NADP/2_6PGL/cl111'
        prmtop_file = 'G6PD_WT_NADP_6PGL_docked.prmtop'
        cpp_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']

        prep_get_frame_from_traj(output_dir, md_traj_dir, prmtop_file, cpp_names)
        res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/run_cpptraj_mmpbsa.sh')

        self.assertEqual(true_res.strip(), res.strip())

    # IMPORTANT!!!
    #  for this test to work you need to comment the line that executes tleap in prepare_mmpbsa_files,
    #  unless tleap is installed in your system
    def test_prep_complex_pdbs_complex_solv(self):
        true_res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/true_res_leaprc_complex_solv')

        mmpbsa_base_dir = 'test_files/test_mmpbsa/prepare_mmpbsa'
        md_traj_dir = '/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa'])
        cl = 'cl111'
        mmpbsa_dir = ''.join([mmpbsa_base_dir, '/', cl])
        pdb_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']
        ligand_params = 'test_files/test_mmpbsa/prepare_mmpbsa/1_G6P'
        leaprc_base = '/'.join([self.base_dir,  'test_files/test_mmpbsa/prepare_mmpbsa/base_leaprc'])
        ligand = 'G6P'

        prep_complex_pdbs(mmpbsa_base_dir, mmpbsa_dir, md_traj_dir, cl, pdb_names, leaprc_base, ligand_params, ligand)
        res = _get_file_content('/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/cl111/0_complex_prep/leaprc_complex_solv']))

        self.assertEqual(true_res.strip(), res.strip())

    # IMPORTANT!!!
    #  for this test to work you need to comment the line that executes tleap in prepare_mmpbsa_files
    #  unless tleap is installed in your system
    def test_prep_complex_pdbs_complex_vac(self):
        true_res = _get_file_content('/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/true_res_leaprc_complex_vac']))

        mmpbsa_base_dir = 'test_files/test_mmpbsa/prepare_mmpbsa'
        md_traj_dir = '/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa'])
        cl = 'cl111'
        mmpbsa_dir = ''.join([mmpbsa_base_dir, '/', cl])
        pdb_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']
        ligand_params = 'test_files/test_mmpbsa/prepare_mmpbsa/1_G6P'
        leaprc_base = '/'.join([self.base_dir,  'test_files/test_mmpbsa/prepare_mmpbsa/base_leaprc'])
        ligand = 'G6P'

        prep_complex_pdbs(mmpbsa_base_dir, mmpbsa_dir, md_traj_dir, cl, pdb_names, leaprc_base, ligand_params, ligand)
        res = _get_file_content('/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/cl111/0_complex_prep/leaprc_complex_vac']))

        self.assertEqual(true_res.strip(), res.strip())

    # IMPORTANT!!!
    #  for this test to work you need to comment the line that executes tleap in prepare_mmpbsa_files
    #  unless tleap is installed in your system
    def test_prep_complex_pdbs_ligand_vac(self):
        true_res = _get_file_content('/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/true_res_leaprc_ligand_vac']))

        mmpbsa_base_dir = 'test_files/test_mmpbsa/process_traj'
        md_traj_dir = '/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa'])
        cl = 'cl111'
        mmpbsa_dir = ''.join([mmpbsa_base_dir, '/', cl])
        pdb_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']
        ligand_params = 'test_files/test_mmpbsa/prepare_mmpbsa/1_G6P'
        leaprc_base = '/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/base_leaprc'])
        ligand = 'G6P'

        prep_complex_pdbs(mmpbsa_base_dir, mmpbsa_dir, md_traj_dir, cl, pdb_names, leaprc_base, ligand_params, ligand)
        res = _get_file_content('/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/cl111/0_complex_prep/leaprc_ligand_vac']))

        self.assertEqual(true_res.strip(), res.strip())

    # IMPORTANT!!!
    #  for this test to work you need to comment the line that executes tleap in prepare_mmpbsa_files
    #  unless tleap is installed in your system
    def test_prep_complex_pdbs_receptor_vac(self):
        true_res = _get_file_content('/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/true_res_leaprc_receptor_vac']))

        mmpbsa_base_dir = 'test_files/test_mmpbsa/process_traj'
        md_traj_dir = '/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa'])
        cl = 'cl111'
        mmpbsa_dir = ''.join([mmpbsa_base_dir, '/', cl])
        pdb_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']
        ligand_params = 'test_files/test_mmpbsa/prepare_mmpbsa/1_G6P'
        leaprc_base = '/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/base_leaprc'])
        ligand = 'G6P'

        prep_complex_pdbs(mmpbsa_base_dir, mmpbsa_dir, md_traj_dir, cl, pdb_names, leaprc_base, ligand_params, ligand)
        res = _get_file_content('/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/cl111/0_complex_prep/leaprc_receptor_vac']))

        self.assertEqual(true_res.strip(), res.strip())


    def test_prep_mmpbsa_in(self):
        true_res = _get_file_content('/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/true_res_mmpbsa_pb.in']))

        frame = 10
        mmpbsa_dict = pb_perframe
        mmpbsa_dict['filename'] = '/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/mmpbsa_pb.in'])
        mmpbsa_dict['start_frame'] = frame
        mmpbsa_dict['end_frame'] = frame

        template_folder = '../../files/jinja_templates'
        mmpbsa_per_frame_template = 'mmpbsa_perframe.jinja'
        prep_mmpbsa_in(template_folder, mmpbsa_per_frame_template, mmpbsa_dict)

        res = _get_file_content('/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/mmpbsa_pb.in']))
        self.assertEqual(true_res.strip(), res.strip())


    def test_prep_mmpbsa_run(self):
        true_res = _get_file_content('test_files/test_mmpbsa/prepare_mmpbsa/true_res_run_mmpbsa.sh')

        mmpbsa_dir = 'test_files/test_mmpbsa/process_traj/cl111'
        mmpbsa_in_file = 'test_files/test_mmpbsa/prepare_mmpbsa/mmpbsa_pb.in'
        md_traj = '/home/marta/G6PD/3_MD_post_dock/1_APO/1_G6P/cl111/G6PD_WT_APO_G6P_docked_cl111_100.crd'

        prep_mmpbsa_run(mmpbsa_dir, mmpbsa_in_file, md_traj)
        res = _get_file_content('/'.join([self.base_dir, 'test_files/test_mmpbsa/prepare_mmpbsa/cl111/run_mmpbsa.sh']))

        self.assertEqual(true_res.strip(), res.strip())


class TestHarvestMMPBSAResults(unittest.TestCase):

    def setUp(self):
        self.base_dir = get_tests_dir()
        os.chdir(self.base_dir)

    def test_harvest_mmpbsa_perframe(self):
        true_res = _get_file_content('test_files/test_mmpbsa/harvest_mmpbsa_results/true_res_mmpbsa_perframe.dat')

        base_dir = 'test_files/test_mmpbsa/harvest_mmpbsa_results'
        n_frames = 9
        file_out = 'per_frame_data'
        harvest_mmpbsa_perframe(base_dir, n_frames, file_out)

        res = _get_file_content('test_files/test_mmpbsa/harvest_mmpbsa_results/mmpbsa_perframe.dat')
        self.assertEqual(true_res.strip(), res.strip())


    def test_average_dG(self):
        true_res = _get_file_content('test_files/test_mmpbsa/harvest_mmpbsa_results/1_2PG/true_res_dG_sumup.dat')

        base_dir = 'test_files/test_mmpbsa/harvest_mmpbsa_results/1_2PG'
        clusters = ['cl000', 'cl100', 'cl110']
        file_out = 'dG_sumup'
        average_dG(base_dir, clusters, file_out)

        res = _get_file_content('test_files/test_mmpbsa/harvest_mmpbsa_results/1_2PG/dG_sumup.dat')
        self.assertEqual(true_res.strip(), res.strip())


    def test_get_ddG(self):
        true_res = _get_file_content('test_files/test_mmpbsa/harvest_mmpbsa_results/true_res_ddG.dat')

        base_dir = 'test_files/test_mmpbsa/harvest_mmpbsa_results/'
        substrate_dir = '1_2PG'
        product_dir = '2_PEP'
        clusters = ['cl000', 'cl100', 'cl110']
        file_out = 'ddG'
        get_ddG(base_dir, substrate_dir, product_dir, clusters, file_out)

        res = _get_file_content('test_files/test_mmpbsa/harvest_mmpbsa_results/ddG.dat')
        self.assertEqual(true_res.strip(), res.strip())

