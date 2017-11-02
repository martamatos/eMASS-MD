def add_active_site_rep(file_vmd, active_site_resids, ligand):
    """
    Given a vmd file, adds the representation for ligand active site residues

    :param file_vmd: path+name for vmd file
    :param active_site_resids: list with active site residues (numbers in strings)
    :param ligand: ligand name
    :return: None
    """

    with open(file_vmd, 'a') as f_out:
        # temp

        # set protein representations
        f_out.write('mol addrep 0\n')
        f_out.write('mol modselect 1 0 protein\n')
        f_out.write('mol modstyle 1 0 NewCartoon 0.300000 10.000000 4.100000 0\n')
        f_out.write('mol modcolor Name\n')

        # set ligand representations
        f_out.write('mol addrep 0\n')
        f_out.write(''.join(['mol modselect 2 0 resname ', ligand, '\n']))
        f_out.write('mol modstyle 2 0 CPK 1.000000 0.300000 12.000000 12.000000\n')
        f_out.write('mol modcolor Name \n')




        f_out.write('mol addrep 0\n')
        f_out.write(''.join(['mol modselect 3 0 resid ', ' '.join(active_site_resids), '\n']))
        f_out.write('mol modstyle 3 0 CPK 1.000000 0.300000 12.000000 12.000000\n')
        f_out.write('mol modcolor 3 0 ResID\n')

        # hide original representation - temp
        f_out.write('mol showrep 0 0 0\n')


def add_vdw_rep_around_ligand(file_vmd, ligand):
    return 0


enzyme = 'GAPDH'
form = 'NAD_remG3P'
form_folder = '_'.join(['3', form])
ligand = '13DPG'
ligand_folder = '_'.join(['2', ligand])

base_form = '_'.join([enzyme, 'WT', form,  ligand, 'docked'])
base_folder = '/'.join(['/home/mrama/Desktop/MD', enzyme, '3_MD_post_dock', form_folder, ligand_folder])

cl_list = ['cl100', 'cl101', 'cl111']

for cl in cl_list:
    to_folder = '/'.join([base_folder, cl])
    file_vmd = '/'.join([to_folder, 'load_traj.tcl'])
    #G6PD
    #active_site_resids = ['235', '177', '147', '238', '181', '239', '178', '234', '377', '346', '339', '222', '215', '344']
    #GAPDH
    active_site_resids = ['149', '148', '150', '208', '209', '231', '176']

    #add_active_site_rep(file_vmd, active_site_resids, ligand)
    add_active_site_rep(file_vmd, active_site_resids, 'DPG')
