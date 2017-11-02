import os
from paramiko import SSHClient

from scp import SCPClient

from src.process_traj.analyze_traj import process_out_files, process_out_files_min, dist_backbone_ligand, write_vmd_script
from src.process_traj.rmsd_MDAnalysis import calculate_traj_rmsd_fit, average_rmsd_from_file
from src.process_traj.traj_utils import convert_traj, merge_trajectories, convert_rst_to_pdb


def progress(filename, size, sent):
    print filename + " " + str(size) + " " + str(sent)


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


def convert_all_rst_to_pdb(base_dir, prmtop):
    for file in os.listdir(base_dir):
        if file.endswith('.rst'):
            convert_rst_to_pdb(file, prmtop)


if __name__ == '__main__':
    enzyme = 'ENO_AB'
    form = 'MG'
    form_folder = '_'.join(['1', form])
    ligand = 'PEP'
    #ligand_list = ['NPD', 'G6P']
    ligand_folder = '_'.join(['2', ligand])
    #ligand_folder = ''
    # ADK1 - AMP
    #binding_residues = ['35', '87', '155', '56']
    # ADK1 - ATP
    #binding_residues = ['122', '12', '9', '14', '13']
    # ENO
    binding_residues = ['340', '391', '369', '207', '370', '368', '367', '338']
    # G6PD1
    #binding_residues = ['234', '176', '146', '237', '180', '238', '177', '233', '376', '345', '338', '221', '214', '343']
    # GAPDH
    #binding_residues = ['148', '147', '149', '207', '208', '230', '175']
    # TALB
    #binding_residues = ['179', '226', '224', '130', '152', '174', '15', '33', '240', '30', '93']
    #binding_residues = None


    md_numbered = '1'  # None or string

    #base_form = '_'.join([enzyme, 'WT', form, '13DPG', 'docked'])
    base_form = '_'.join(['ENO_AB', 'WT', 'APO', ligand, 'docked'])
    #base_form = '_'.join(['TALB_WT_HALO_S7P_remS7P_F6P_docked'])
    base_folder = '/'.join(['/home/mrama/Desktop/MD', enzyme, '3_MD_post_dock2016', form_folder, ligand_folder])
    #base_folder = '/'.join(['/home/mrama/Desktop/MD', enzyme, '3_MD_post_dock2016', form_folder])

    cl_list = ['cl302', 'cl312', 'cl402']
    #server = 'sbrgstructure.dynamic.ucsd.edu'
    server = 'ssb3.ucsd.edu'
    traj_format_list = ['crd', 'crd', 'crd']
    traj_type = '100'

    copy_from_server = False

    #if ligand == '13DPG':
    #    ligand = 'DPG'
    #elif ligand == '6PGL':
    #    ligand = '6PG'
    #ligand = 'S7P'
    for i, cl in enumerate(cl_list):
        traj_format = traj_format_list[i]
        from_folder = '/'.join(['/home/marta', enzyme, '3_MD_post_dock2016', form_folder, ligand_folder, cl])
        #from_folder = '/'.join(['/home/marta', enzyme, '3_MD_post_dock2016', form_folder, cl])
        to_folder = '/'.join([base_folder, cl])


        # download trajectories
        if not os.path.isdir(to_folder):
            os.mkdir(to_folder)
        if copy_from_server:
            scp_traj_files(server, from_folder, to_folder, traj_format)

        traj = ''.join([to_folder, '/', base_form, '_', cl, '_', traj_type, '.', traj_format]) if not md_numbered \
            else ''.join([to_folder, '/',  base_form, '_', cl, '_', md_numbered, '_', traj_type, '.', traj_format])
        prmtop = ''.join([to_folder, '/',  base_form, '.prmtop'])

        #convert trajectory format from crd to dcd
        if traj_format == 'crd':
            convert_traj(traj, prmtop, format_out='dcd', remove_orig=True)
            traj = ''.join([to_folder, '/', base_form, '_', cl, '_', traj_type, '.dcd']) if not md_numbered \
                else ''.join([to_folder, '/',  base_form, '_', cl, '_', md_numbered, '_', traj_type, '.dcd'])
            traj_format = 'dcd'

        # process out files
        '''
        file_in = ''.join([to_folder, '/', base_form, '_min_v0_l1.out'])
        file_out = ''.join([to_folder, '/data/min'])
        process_out_files_min(to_folder, file_in, file_out)
        file_in = ''.join([to_folder, '/', base_form, '_equil_v0_l1.out'])
        file_out = ''.join([to_folder, '/data/equil'])
        process_out_files(to_folder, file_in, file_out)
        file_in = ''.join([to_folder, '/', base_form, '_prod_v0_l1.out'])
        file_out = ''.join([to_folder, '/data/prod'])
        process_out_files(to_folder, file_in, file_out)
        file_in = ''.join([to_folder, '/', base_form, '_prod_v0_l2.out'])
        file_out = ''.join([to_folder, '/data/prod'])
        #process_out_files(to_folder, file_in, file_out)'''

        if md_numbered != '1' and md_numbered != 'combined':
            n_trajs = int(md_numbered)
            traj_list = []
            for traj_i in range(1, n_trajs+1):
                traj_list.append(''.join([to_folder, '/',  base_form, '_', cl, '_', str(traj_i), '_', traj_type, '.', traj_format]))
            print traj_list
            traj = ''.join([to_folder, '/',  base_form, '_', cl, '_combined_', traj_type, '.', traj_format])
            merge_trajectories(traj_list, prmtop, traj)

        convert_all_rst_to_pdb(to_folder, prmtop)


        # get rmsd
        file_rmsd = ''.join([to_folder, '/rmsd_backbone_', traj_type, '_', cl])
        calculate_traj_rmsd_fit(traj, prmtop, file_rmsd, selection='backbone')
        file_out = ''.join([to_folder, '/rmsd_backbone_average_', traj_type])
        average_rmsd_from_file('.'.join([file_rmsd, 'dat']), file_out)

        ligand_list = ['PEP']
        for ligand in ligand_list:
            file_out = ''.join([to_folder, '/rmsd_', ligand, '_', cl])
            calculate_traj_rmsd_fit(traj, prmtop, file_out, selection=' '.join(['resname', ligand]))

        # get ligand distances
        for ligand in ligand_list:
            file_out = ''.join([to_folder, '/', ligand, '_dist_', cl])
            atom_pairs = dist_backbone_ligand(traj, prmtop, ligand, file_out)

        # write vmd scripts
        ligand = 'PEP'
        file_out = '/'.join([to_folder, 'load_traj.tcl'])
        write_vmd_script(traj, prmtop, file_out, ligand, atom_pairs=atom_pairs, binding_residues=binding_residues)
