import os
from paramiko import SSHClient
from scp import SCPClient
from src.process_traj.analyze_traj import  get_atoms_dist, process_out_files, process_out_files_min, write_vmd_script
from src.process_traj.rmsd_MDAnalysis import calculate_traj_rmsd_fit, average_rmsd_from_file
from src.process_traj.traj_utils import convert_traj, merge_trajectories, recenter_traj


def scp_traj_files(server, from_folder, to_folder, traj_format):
    """
        Copy out, crd, prmtop files from server to local machine

    :param server: server name, e.g. ssb3.ucsd.edu
    :param from_folder: where the files are located
    :param to_folder: where the files should be copied to
    :return: trajctory format extension, either .crd or .dcd
    """

    ssh = SSHClient()
    ssh.load_system_host_keys()
    ssh.connect(server, username='marta')

    with SCPClient(ssh.get_transport(), sanitize=lambda x: x) as scp:
        print 'out files'
        scp.get(''.join([from_folder, '/*.out']), to_folder, preserve_times=True, recursive=True)
        print 'rst files'
        scp.get(''.join([from_folder, '/*.rst']), to_folder, preserve_times=True, recursive=True)
        print 'prmtop file'
        scp.get(''.join([from_folder, '/*.prmtop']), to_folder, preserve_times=True, recursive=True)
        print 'traj file'
        scp.get(''.join([from_folder, '/*100.*']), to_folder, preserve_times=True, recursive=True)


def analyze_ENO_1E9I_APO():

    base_folder = '/home/mrama/Desktop/MD/ENO_1E9I/1_MD/1_APO'
    to_folder = base_folder
    prmtop = '/'.join([base_folder, '1E9I_D_WT_APO.prmtop'])
    traj = '/'.join([base_folder, '1E9I_D_WT_APO_100_l3.crd'])

    # process out files
    equil_time = None
    prod_time = None
    base_name = '1E9I_D_WT_APO'
    #process_out_files(base_folder, base_name, equil_time, prod_time)

    rst_file_out = '1E9I_D_WT_APO_prod_v0_l2_recentered.rst'
    n_residues = 433
    #recenter_traj(traj, prmtop, rst_file_out, n_residues)

    traj_list = ['1E9I_D_WT_APO_100_l1.dcd', '1E9I_D_WT_APO_100_l2.dcd', '1E9I_D_WT_APO_100_l3.dcd']
    traj_list = ['/'.join([base_folder, traj]) for traj in traj_list]
    traj_out = '/'.join([base_folder, '1E9I_D_WT_APO_100_merged.dcd'])
    #merge_trajectories(traj_list, prmtop, traj_out)
    traj = traj_out


    # process out files
    file_in = ''.join([to_folder, '/', base_name, '_min_v0_l1.out'])
    file_out = ''.join([to_folder, '/data/min'])
    process_out_files_min(to_folder, file_in, file_out)
    file_in = ''.join([to_folder, '/', base_name, '_equil_v0_l1.out'])
    file_out = ''.join([to_folder, '/data/equil'])
    process_out_files(to_folder, file_in, file_out)
    file_in = ''.join([to_folder, '/', base_name, '_prod_v0_l2.out'])
    file_out = ''.join([to_folder, '/data/prod'])
    process_out_files(to_folder, file_in, file_out)
    file_in = ''.join([to_folder, '/', base_name, '_prod_v0_l3.out'])
    file_out = ''.join([to_folder, '/data/prod'])
    process_out_files(to_folder, file_in, file_out)

    # convert trajectory
    #convert_traj(traj, prmtop, format_out='dcd', remove_orig=True, use_cpptraj=True)
    #traj = '/'.join([base_folder, '1E9I_D_WT_APO_100_l3.dcd'])

    # get rmsd
    file_rmsd = ''.join([base_folder, '/rmsd_backbone_merged'])
    calculate_traj_rmsd_fit(traj, prmtop, file_rmsd, selection='backbone')
    file_out = ''.join([base_folder, '/rmsd_backbone_average_merged'])
    average_rmsd_from_file('.'.join([file_rmsd, 'dat']), file_out)

    # get atoms distance
    atom_pairs = [(6418, 3573), (6418, 4257), (6418, 4658), (6419, 4658), (6419, 5064), (6419, 5505), (6419, 3573)]
    file_out = ''.join([base_folder, '/MG_restraints_merged'])
    get_atoms_dist(traj, prmtop, atom_pairs, file_out)

    # write vmd scripts
    file_out = '/'.join([base_folder, 'load_traj_merged.tcl'])
    write_vmd_script(traj, prmtop, file_out)



def analyze_GAPDH_APO():

    base_folder = '/home/mrama/Desktop/MD/GAPDH/1_MD/1_APO'
    prmtop = '/'.join([base_folder, 'GAPDH_WT_APO.prmtop'])
    traj = '/'.join([base_folder, 'GAPDH_WT_APO_100.dcd'])


    # process out files
    equil_time = None
    prod_time = None
    base_name = 'TALB_WT_HALO'
    #process_out_files(base_folder, base_name, equil_time, prod_time)

    # get rmsd
    file_rmsd = ''.join([base_folder, '/rmsd_backbone'])
    calculate_traj_rmsd_fit(traj, prmtop, file_rmsd, selection='backbone')
    file_out = ''.join([base_folder, '/rmsd_backbone_average'])
    average_rmsd_from_file('.'.join([file_rmsd, 'dat']), file_out)

    # write vmd scripts
    file_out = '/'.join([base_folder, 'load_traj.tcl'])
    write_vmd_script(traj, prmtop, file_out)


def analyze_GAPDH_NAD():

    base_folder = '/home/mrama/Desktop/MD/GAPDH/1_MD/2_NAD'
    prmtop = '/'.join([base_folder, 'GAPDH_WT_NAD.prmtop'])
    traj = '/'.join([base_folder, 'GAPDH_WT_NAD_100.dcd'])


    # process out files
    equil_time = None
    prod_time = None
    base_name = 'TALB_WT_HALO'
    #process_out_files(base_folder, base_name, equil_time, prod_time)

    # get rmsd
    file_rmsd = ''.join([base_folder, '/rmsd_backbone'])
    calculate_traj_rmsd_fit(traj, prmtop, file_rmsd, selection='backbone')
    file_out = ''.join([base_folder, '/rmsd_backbone_average'])
    average_rmsd_from_file('.'.join([file_rmsd, 'dat']), file_out)

    # write vmd scripts
    file_out = '/'.join([base_folder, 'load_traj.tcl'])
    write_vmd_script(traj, prmtop, file_out)


def analyze_GAPDH_HALO_G3P():

    base_folder = '/home/mrama/Desktop/MD/GAPDH/1_MD/3_HALO'
    prmtop = '/'.join([base_folder, 'GAPDH_WT_HALO_G3P.prmtop'])
    traj = '/'.join([base_folder, 'GAPDH_WT_HALO_G3P_100.dcd'])


    # process out files
    equil_time = None
    prod_time = None
    base_name = 'TALB_WT_HALO'
    #process_out_files(base_folder, base_name, equil_time, prod_time)

    # get rmsd
    file_rmsd = ''.join([base_folder, '/rmsd_backbone'])
    calculate_traj_rmsd_fit(traj, prmtop, file_rmsd, selection='backbone')
    file_out = ''.join([base_folder, '/rmsd_backbone_average'])
    average_rmsd_from_file('.'.join([file_rmsd, 'dat']), file_out)

    # write vmd scripts
    file_out = '/'.join([base_folder, 'load_traj.tcl'])
    write_vmd_script(traj, prmtop, file_out)


def analyze_TALB_APO():

    base_folder = '/home/mrama/Desktop/MD/TALB/1_MD/1_APO'
    prmtop = '/'.join([base_folder, '1UCW_WT.prmtop'])
    traj = '/'.join([base_folder, 'TALB_WT_APO_100.dcd'])


    # process out files
    equil_time = None
    prod_time = None
    base_name = 'TALB_WT_APO'
    #process_out_files(base_folder, base_name, equil_time, prod_time)

    # get rmsd
    file_rmsd = ''.join([base_folder, '/rmsd_backbone'])
    calculate_traj_rmsd_fit(traj, prmtop, file_rmsd, selection='backbone')
    file_out = ''.join([base_folder, '/rmsd_backbone_average'])
    average_rmsd_from_file('.'.join([file_rmsd, 'dat']), file_out)

    # write vmd scripts
    file_out = '/'.join([base_folder, 'load_traj.tcl'])
    write_vmd_script(traj, prmtop, file_out)


def analyze_TALB_HALO():

    base_folder = '/home/mrama/Desktop/MD/TALB/1_MD/2_HALO'
    prmtop = '/'.join([base_folder, '1UCW_WT_S7P.prmtop'])
    traj = '/'.join([base_folder, 'TALB_WT_HALO_S7P_100.dcd'])


    # process out files
    equil_time = None
    prod_time = None
    base_name = 'TALB_WT_HALO'
    #process_out_files(base_folder, base_name, equil_time, prod_time)

    # get rmsd
    file_rmsd = ''.join([base_folder, '/rmsd_backbone'])
    calculate_traj_rmsd_fit(traj, prmtop, file_rmsd, selection='backbone')
    file_out = ''.join([base_folder, '/rmsd_backbone_average'])
    average_rmsd_from_file('.'.join([file_rmsd, 'dat']), file_out)

    # write vmd scripts
    file_out = '/'.join([base_folder, 'load_traj.tcl'])
    write_vmd_script(traj, prmtop, file_out)


if __name__ == '__main__':
    #analyze_TALB_APO()
    #analyze_TALB_HALO()
    #analyze_GAPDH_APO()
    #analyze_GAPDH_NAD()
    #analyze_GAPDH_HALO_G3P()
    #analyze_ENO_1E9I_APO()
    analyze_ENO_1E9I_APO()