import os
import subprocess

import mdtraj


CPPTRAJ_DIR = '/home/mrama/Programs/amber14'


def weed_traj(md_traj_dir, script_dir, base_name, frame_intervals, md_sim, cl):
    """
        Creates the cpptraj files needed to weed an MD trajectory and runs them.
        In the end the mdcrd file is deleted.

    :param md_traj_dir: absolute path to the folder where the md simulations are, ignoring the cluster,
                        e.g.: /home/marta/G6PD/3_MD_post_dock/1_APO/1_G6P
    :param script_dir: absolute path to where the cpptraj script is located
    :param base_name: base name of the MD simulation file
    :param frame_intervals: list with the intervals between each frame that should be saved
    :param md_sim: in case of more than one run, which one should be chosen
    :param cl: name of which cluster is the MD simulation from, eg. cl100
    :return: None
    # """

    for frame_interval in frame_intervals:
        # created cpptraj script
        with open('/'.join([script_dir, 'weed_traj.in']), 'w') as f_weed_traj:
            f_weed_traj.write(''.join(['trajin mdcrd 1 last ', str(frame_interval), '\n']))
            f_weed_traj.write(''.join(['trajout ', base_name, '_', cl, '_', str(md_sim), '_', str(frame_interval), '.crd\n']))

        # go to the folder where the md simulation is
        os.chdir(md_traj_dir)
        # run cpptraj script
        subprocess.call(''.join(['$AMBERHOME/bin/cpptraj -i ', script_dir, '/weed_traj.in -p ', base_name]), shell=True)
        subprocess.call('rm -rf mdcrd', shell=True)


def convert_traj(traj, prmtop, format_out, remove_orig=False, use_cpptraj=False):
    """
        Convert an MD trajectory to a given format, uses mdtraj or cpptraj, depending on the value of use_cpptraj.

    :param traj: path to MD trajectory
    :param prmtop: path to prmtop file
    :param format_out: string with the output format, e.g. 'dcd'
    :param remove_orig: a boolean specifying whether or not the original trajectory should be deleted
    :param use_cpptraj: a boolean specifying whether mdtraj or cpptraj should be used to convert the trajectory format
    :return: None
    """
    if not use_cpptraj:
        t = mdtraj.load(traj, top=prmtop)
        t.save_dcd(''.join([traj[:-3], format_out]))
    else:
        traj_dir = traj.split('/')
        traj_dir = '/'. join(traj_dir[:-1])

        with open('/'.join([traj_dir, 'weed_traj.in']), 'w') as f_weed_traj:
            f_weed_traj.write(''.join(['trajin ', traj, '\n']))
            f_weed_traj.write(''.join(['trajout ', traj[:-3], format_out, '\n']))

        # run cpptraj script
        subprocess.call(''.join(['$AMBERHOME/bin/cpptraj -i ', traj_dir, '/weed_traj.in -p ', prmtop]), shell=True)

    if remove_orig and os.path.isfile(''.join([traj[:-3], format_out])):
        os.remove(traj)


def merge_trajectories(traj_list, prmtop, traj_out):
    """
        Given a list of MD trajectories, concatenates them all into a single one.

    :param traj_list: a list of paths to MD trajectories
    :param prmtop: path to prmtop file
    :param traj_out: path to merged trajectory
    :return: None
    """
    t = mdtraj.load(traj_list[0], top=prmtop)
    for traj_i in range(1, len(traj_list)):
        temp_traj = mdtraj.load(traj_list[traj_i], top=prmtop)
        t = t.join(temp_traj)
    t.save_dcd(traj_out)


def recenter_traj(traj, prmtop, rst_file_out, n_residues):
    """
        Given a trajectory and prmtop file, generates a restard file where the system is re-centered.

    :param traj: path to the MD trajectory
    :param prmtop: path to the prmtop file
    :param rst_file_out: path to the resulting restart file
    :param n_residues: number of residues in the system
    :return: None
    """

    traj_temp = traj.split('/')
    traj_dir = '/'. join(traj_temp[:-1])
    traj_name = traj_temp[-1:][0]


    with open('/'.join([traj_dir, 'recenter_traj.in']), 'w') as f_recenter:
        f_recenter.write(''.join(['trajin ', traj_name, ' lastframe\n']))
        f_recenter.write(''.join(['trajout ', rst_file_out, ' restart\n']))
        f_recenter.write(''.join(['center :1-', str(n_residues), '\n']))
        f_recenter.write('image\n')
        f_recenter.write('go')

    subprocess.call(''.join(['$AMBERHOME/bin/cpptraj -i ', traj_dir, '/recenter_traj.in -p ', prmtop]), shell=True)


def convert_rst_to_pdb(rst_file, prmtop, cpp_traj_file=None):
    """
    Given a restart file and the respective prmtop, converts it into a pdb file.

    :param rst_file: path to the rst file to be converted
    :param prmtop: path to the prmtop file
    :param cpp_traj_file: path to the cpptraj file to be produced, optional.
    :return: None
    """

    file_out = cpp_traj_file if cpp_traj_file else './temp.in'
    with open(file_out, 'w') as f_out:
        f_out.write(''.join(['trajin ', rst_file, ' lastframe\n']))
        f_out.write(''.join(['trajout ', rst_file[:-3], 'pdb pdb\n']))
        f_out.write('go\n')

    subprocess.call(''.join(['$AMBERHOME/bin/cpptraj -i ', file_out, ' -p ', prmtop]), shell=True)
    if not cpp_traj_file:
        os.remove(file_out)


def extract_pdb_frames_from_traj(prmtop, traj_file, frame_list, base_file_out, ligand):
    """
    Given a prmtop+trajectory path and a list frames, extracts this frames to pdb files with name
    base_file_out_framenumber.

    :param prmtop: path+name for prmtop file
    :param traj_file: MD trajectory file path+name
    :param frame_list: list with frames to be extracted from trajectory
    :param base_file_out: base path+name for output file
    :return: None
    """

    for frame_i in frame_list:
        file_out = './temp.in'
        with open(file_out, 'w') as f_out:
            f_out.write(''.join(['trajin ', traj_file, ' ', str(frame_i), ' ', str(frame_i), '\n']))
            f_out.write(''.join(['trajout ', base_file_out, '-%.3d_receptor_' % (frame_i-249), ligand,  '_docked.pdb pdb\n']))
            f_out.write('strip :WAT:Na+\n')
            f_out.write('go\n')

        subprocess.call(''.join(['$AMBERHOME/bin/cpptraj -i ', file_out, ' -p ', prmtop]), shell=True)
        os.remove(file_out)



prmtop = '/home/mrama/Desktop/MD/ENO_AB/3_MD_post_dock/1_MG/1_2PG/cl000/ENO_AB_WT_MG_2PG_docked.prmtop'
traj_file = '/home/mrama/Desktop/MD/ENO_AB/3_MD_post_dock/1_MG/1_2PG/cl000/ENO_AB_WT_MG_2PG_docked_cl000_merged_100.dcd'
frame_list = range(250, 502)
ligand = '2PG'
base_file_out = '/home/mrama/Desktop/MD/ENO_AB/2_Dock/1_MG/1_2PG/dock_prep/ENO_AB_WT_MG'
extract_pdb_frames_from_traj(prmtop, traj_file, frame_list, base_file_out, ligand)