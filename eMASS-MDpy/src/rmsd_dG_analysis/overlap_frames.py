

def create_vmd_file(file_names_list, file_out, ligand_list, binding_residues):
    """
    Creates .tcl file to be loaded by vmd with all frames in file_names_list.

    :param file_names_list: list with paths to pdb files
    :param file_out: path for .tcl file
    :param ligand_list: a list of ligands bound to the enzyme in file_names_list
    :param binding_residues: a list with the residue numbers
    :return: None
    """

    with open(file_out, 'w') as f_out:
        for frame_i, file_name in enumerate(file_names_list):
            f_out.write(''.join(['mol new {', file_name, '} type {pdb} first 0 last -1 step 1 waitfor 1\n']))



            f_out.write(''.join(['mol addrep ', str(frame_i), '\n']))
            f_out.write(''.join(['mol modselect 1 ', str(frame_i), ' protein\n']))
            f_out.write(''.join(['mol modstyle 1 ', str(frame_i), ' NewCartoon 0.300000 10.000000 4.100000 0\n']))
            f_out.write('mol modcolor Name\n')
            f_out.write('\n')

            for ligand_i, ligand in enumerate(ligand_list):
                f_out.write(''.join(['mol addrep ', str(frame_i), '\n']))
                f_out.write(''.join(['mol modselect ', str(ligand_i+2), ' ', str(frame_i), ' resname ', ligand, '\n']))
                f_out.write(''.join(['mol modstyle ', str(ligand_i+2), ' ', str(frame_i), ' Licorice 0.200000 12.000000 12.000000\n']))

                if ligand_i == 0:
                    f_out.write('mol modcolor Name\n')
                else:
                    f_out.write(''.join(['mol modcolor ', str(ligand_i+2), ' ', str(frame_i), ' ColorID ', str(ligand_i+1), '\n']))
                f_out.write('\n')


            f_out.write(''.join(['mol addrep ', str(frame_i), '\n']))
            f_out.write(''.join(['mol modselect ', str(len(ligand_list)+2), ' ', str(frame_i), ' residue ', ' '.join(binding_residues), '\n']))
            f_out.write(''.join(['mol modstyle ', str(len(ligand_list)+2), ' ', str(frame_i), ' VDW 0.6 12.000000\n']))
            f_out.write(''.join(['mol modcolor ', str(len(ligand_list)+2), ' ', str(frame_i), ' ColorID 11\n']))

            f_out.write(''.join(['mol showrep ', str(frame_i), ' 0 0\n']))
            f_out.write('\n\n')

        f_out.write('color Display Background white\n')


def create_tcl_for_eno_1e9i(base_dir):
    """
    Defines all arguments for create_vmd_file(), to create .tcl files with all starting frames for the MD simulations
    with ENO-2PG and ENO-PEP (ENO structure 1E9I).
    :param base_dir: path to folder with MD simulations
    :return: None
    """

    binding_residues = ['340', '391', '369', '207', '370', '368', '367', '338']
    base_dir_ligand = ''.join([base_dir, '/ENO_1E9I/3_MD_post_dock2016/1_APO/1_2PG/'])
    file_names_list = [''.join([base_dir_ligand, 'cl000/0_complex_prep/1E9I_WT_APO-253_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl000_2/0_complex_prep/1E9I_WT_APO-093_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl000_3/0_complex_prep/1E9I_WT_APO-205_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl000_4/0_complex_prep/1E9I_WT_APO-738_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl100/0_complex_prep/1E9I_WT_APO-094_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl100_2/0_complex_prep/1E9I_WT_APO-744_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl100_3/0_complex_prep/1E9I_WT_APO-384_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl100_4/0_complex_prep/1E9I_WT_APO-504_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl200/0_complex_prep/1E9I_WT_APO-268_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl200_2/0_complex_prep/1E9I_WT_APO-405_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl200_3/0_complex_prep/1E9I_WT_APO-715_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl200_4/0_complex_prep/1E9I_WT_APO-308_receptor_2PG_docked.pdb'])]

    file_out = ''.join([base_dir, '/ENO_1E9I/3_MD_post_dock2016/load_2pg_frames.tcl'])
    ligand_list = ['2PG', 'MG']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)

    base_dir_ligand = ''.join([base_dir, '/ENO_1E9I/3_MD_post_dock2016/1_APO/2_PEP/'])
    file_names_list = [''.join([base_dir_ligand, 'cl000/0_complex_prep/1E9I_WT_APO-241_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl000_2/0_complex_prep/1E9I_WT_APO-040_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl000_3/0_complex_prep/1E9I_WT_APO-567_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl000_4/0_complex_prep/1E9I_WT_APO-138_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl100/0_complex_prep/1E9I_WT_APO-228_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl100_2/0_complex_prep/1E9I_WT_APO-740_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl100_3/0_complex_prep/1E9I_WT_APO-353_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl100_4/0_complex_prep/1E9I_WT_APO-431_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl110/0_complex_prep/1E9I_WT_APO-178_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl110_2/0_complex_prep/1E9I_WT_APO-468_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl110_3/0_complex_prep/1E9I_WT_APO-079_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl110_4/0_complex_prep/1E9I_WT_APO-612_receptor_PEP_docked.pdb'])]
                       #'/home/mrama/Desktop/MD/ENO_1E9I/0_Crystal_structures/3H8A_A_WT_PEP.pdb']

    file_out = ''.join([base_dir, '/ENO_1E9I/3_MD_post_dock2016/load_pep_frames.tcl'])
    ligand_list = ['PEP', 'MG']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)


def create_tcl_for_eno_ab(base_dir):
    """
    Defines all arguments for create_vmd_file(), to create .tcl files with all starting frames for the MD simulations
    with ENO-2PG and ENO-PEP (ENO structure: 3H8A).
    :return: None
    """

    binding_residues = ['340', '391', '369', '207', '370', '368', '367', '338']

    base_dir_ligand = ''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/1_MG/1_2PG/'])
    file_names_list = [''.join([base_dir_ligand, 'cl102/0_complex_prep/ENO_AB_WT_MG-047_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl202/0_complex_prep/ENO_AB_WT_MG-109_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl203/0_complex_prep/ENO_AB_WT_MG-175_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl000/0_complex_prep/ENO_AB_WT_MG-001_receptor_2PG_docked.pdb']),
                       ''.join([base_dir, '/ENO_1E9I/0_Crystal_structures/3H8A_A_WT_2PG.pdb'])]

    file_out = ''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/load_2pg_frames.tcl'])
    ligand_list = ['2PG', 'MG']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)

    base_dir_ligand = ''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/1_MG/2_PEP/'])
    file_names_list = [''.join([base_dir_ligand, 'cl302/0_complex_prep/ENO_AB_WT_MG-121_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl312/0_complex_prep/ENO_AB_WT_MG-244_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl402/0_complex_prep/ENO_AB_WT_MG-058_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl000/0_complex_prep/ENO_AB_WT_MG-001_receptor_PEP_docked.pdb']),
                       ''.join([base_dir, '/ENO_1E9I/0_Crystal_structures/3H8A_A_WT_PEP.pdb'])]

    file_out = ''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/load_pep_frames.tcl'])
    ligand_list = ['PEP', 'MG']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)


def create_tcl_for_gapd(base_dir):
    """
    Defines all arguments for create_vmd_file(), to create .tcl files with all starting frames for the MD simulations
    with GAPD-NAD-G3P and GAPD-NADH-DPG.
    :return: None
    """

    binding_residues = ['148', '147', '149', '207', '208', '230', '175']

    base_dir_ligand = ''.join([base_dir, '/GAPDH/3_MD_post_dock2016/'])
    file_names_list = [''.join([base_dir_ligand, '1_NAD/1_G3P/cl000/0_complex_prep/GAPDH_WT_NAD-076_receptor_G3P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/1_G3P/cl000_2/0_complex_prep/GAPDH_WT_NAD-252_receptor_G3P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/1_G3P/cl001/0_complex_prep/GAPDH_WT_NAD-068_receptor_G3P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/1_G3P/cl001_2/0_complex_prep/GAPDH_WT_NAD-343_receptor_G3P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/1_G3P/cl012/0_complex_prep/GAPDH_WT_NAD-159_receptor_G3P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/1_G3P/cl012_2/0_complex_prep/GAPDH_WT_NAD-540_receptor_G3P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl010/0_complex_prep/GAPDH_WT_NAD_rem3PG-275_receptor_G3P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl010_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-301_receptor_G3P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl100/0_complex_prep/GAPDH_WT_NAD_rem3PG-626_receptor_G3P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl100_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-107_receptor_G3P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200/0_complex_prep/GAPDH_WT_NAD_rem3PG-025_receptor_G3P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-430_receptor_G3P_docked.pdb']),
                       ''.join([base_dir, '/GAPDH/0_Crystal_structures/1DC4_chainA.pdb'])]

    file_out = ''.join([base_dir, '/GAPDH/3_MD_post_dock2016/load_g3p_frames.tcl'])
    ligand_list = ['G3P', 'NAD']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)


    base_dir_ligand =''.join([base_dir, '/GAPDH/3_MD_post_dock2016/'])
    file_names_list = [''.join([base_dir_ligand, '1_NAD/2_DPG/cl100/0_complex_prep/GAPDH_WT_NAD-027_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl100_2/0_complex_prep/GAPDH_WT_NAD-139_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl101/0_complex_prep/GAPDH_WT_NAD-592_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl101_2/0_complex_prep/GAPDH_WT_NAD-234_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl220/0_complex_prep/GAPDH_WT_NAD-316_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl220_2/0_complex_prep/GAPDH_WT_NAD-399_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl100/0_complex_prep/GAPDH_WT_NAD_rem3PG-010_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl100_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-640_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl120/0_complex_prep/GAPDH_WT_NAD_rem3PG-242_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl120_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-348_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl210/0_complex_prep/GAPDH_WT_NAD_rem3PG-516_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl210_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-140_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir, '/GAPDH/0_Crystal_structures/1DC4_chainA.pdb'])]


    file_out = ''.join([base_dir, '/GAPDH/3_MD_post_dock2016/load_dpg_frames.tcl'])
    ligand_list = ['DPG', 'NAD']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)


def create_tcl_for_talb(base_dir):
    """
    Defines all arguments for create_vmd_file(), to create .tcl files with all starting frames for the MD simulations
    with TALB-S7P and TALB-F6P.
    :return: None
    """

    binding_residues = ['179', '226', '224', '130', '152', '174', '15', '33', '240', '30', '93']
    base_dir_ligand = ''.join([base_dir, '/TalB/3_MD_post_dock2016/'])
    file_names_list = [''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl101/0_complex_prep/TALB_WT_APO-103_receptor_S7P_linear_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl101_2/0_complex_prep/TALB_WT_APO-515_receptor_S7P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl102/0_complex_prep/TALB_WT_APO-091_receptor_S7P_linear_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl102_2/0_complex_prep/TALB_WT_APO-535_receptor_S7P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl212/0_complex_prep/TALB_WT_APO-036_receptor_S7P_linear_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl212_2/0_complex_prep/TALB_WT_APO-777_receptor_S7P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl101/0_complex_prep/TALB_WT_HALO_S7P_remS7P-237_receptor_S7P_linear_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl101_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-260_receptor_S7P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl302/0_complex_prep/TALB_WT_HALO_S7P_remS7P-091_receptor_S7P_linear_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl302_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-492_receptor_S7P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl322/0_complex_prep/TALB_WT_HALO_S7P_remS7P-042_receptor_S7P_linear_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl322_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-184_receptor_S7P_docked.pdb']),
                       ''.join([base_dir, '/TalB/0_Crystalline_structures/1UCW_wt_S7P_watbox.pdb'])]

    file_out = ''.join([base_dir, '/TalB/3_MD_post_dock2016/load_s7p_frames.tcl'])
    ligand_list = ['S7P']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)


    base_dir_ligand = ''.join([base_dir, '/TalB/3_MD_post_dock2016/'])
    file_names_list = [''.join([base_dir_ligand, '1_APO/2_F6P/cl202/0_complex_prep/TALB_WT_APO-533_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl202_2/0_complex_prep/TALB_WT_APO-314_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl212/0_complex_prep/TALB_WT_APO-139_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl212_2/0_complex_prep/TALB_WT_APO-464_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl321/0_complex_prep/TALB_WT_APO-799_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl321_2/0_complex_prep/TALB_WT_APO-003_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl101/0_complex_prep/TALB_WT_HALO_S7P_remS7P-492_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl101_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-056_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl201/0_complex_prep/TALB_WT_HALO_S7P_remS7P-381_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl201_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-117_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl311/0_complex_prep/TALB_WT_HALO_S7P_remS7P-014_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl311_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-222_receptor_F6P_docked.pdb']),
                       ''.join([base_dir, '/TalB/0_Crystalline_structures/1UCW_wt_S7P_watbox.pdb'])]


    file_out = ''.join([base_dir, '/TalB/3_MD_post_dock2016/load_f6p_frames.tcl'])
    ligand_list = ['F6P']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)


def main():
    base_dir = '/home/mrama/Dropbox/PhD_stuff/Projects/MD/eMASS-MD/MD_data/'
    create_tcl_for_eno_1e9i(base_dir)
    create_tcl_for_eno_ab(base_dir)
    create_tcl_for_gapd(base_dir)
    create_tcl_for_talb(base_dir)

if __name__ == '__main__':
    main()

