

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
                       ''.join([base_dir_ligand, 'cl102_2/0_complex_prep/ENO_AB_WT_MG-004_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl102_3/0_complex_prep/ENO_AB_WT_MG-059_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl102_4/0_complex_prep/ENO_AB_WT_MG-064_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl102_5/0_complex_prep/ENO_AB_WT_MG-091_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl102_6/0_complex_prep/ENO_AB_WT_MG-053_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl102_7/0_complex_prep/ENO_AB_WT_MG-054_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl102_8/0_complex_prep/ENO_AB_WT_MG-073_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl102_9/0_complex_prep/ENO_AB_WT_MG-078_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl102_10/0_complex_prep/ENO_AB_WT_MG-090_receptor_2PG_docked.pdb']),

                       ''.join([base_dir_ligand, 'cl202/0_complex_prep/ENO_AB_WT_MG-109_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl202_2/0_complex_prep/ENO_AB_WT_MG-015_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl202_3/0_complex_prep/ENO_AB_WT_MG-038_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl202_4/0_complex_prep/ENO_AB_WT_MG-130_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl202_5/0_complex_prep/ENO_AB_WT_MG-153_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl202_6/0_complex_prep/ENO_AB_WT_MG-029_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl202_7/0_complex_prep/ENO_AB_WT_MG-096_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl202_8/0_complex_prep/ENO_AB_WT_MG-105_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl202_9/0_complex_prep/ENO_AB_WT_MG-123_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl202_10/0_complex_prep/ENO_AB_WT_MG-144_receptor_2PG_docked.pdb']),

                       ''.join([base_dir_ligand, 'cl203/0_complex_prep/ENO_AB_WT_MG-175_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl203_2/0_complex_prep/ENO_AB_WT_MG-003_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl203_3/0_complex_prep/ENO_AB_WT_MG-148_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl203_4/0_complex_prep/ENO_AB_WT_MG-159_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl203_5/0_complex_prep/ENO_AB_WT_MG-161_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl203_6/0_complex_prep/ENO_AB_WT_MG-017_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl203_7/0_complex_prep/ENO_AB_WT_MG-127_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl203_8/0_complex_prep/ENO_AB_WT_MG-179_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl203_9/0_complex_prep/ENO_AB_WT_MG-103_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl203_10/0_complex_prep/ENO_AB_WT_MG-167_receptor_2PG_docked.pdb']),

                       ''.join([base_dir_ligand, 'cl103/0_complex_prep/ENO_AB_WT_MG-001_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl103_2/0_complex_prep/ENO_AB_WT_MG-048_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl103_3/0_complex_prep/ENO_AB_WT_MG-060_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl103_4/0_complex_prep/ENO_AB_WT_MG-094_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl103_5/0_complex_prep/ENO_AB_WT_MG-162_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl103_6/0_complex_prep/ENO_AB_WT_MG-050_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl103_7/0_complex_prep/ENO_AB_WT_MG-057_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl103_8/0_complex_prep/ENO_AB_WT_MG-079_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl103_9/0_complex_prep/ENO_AB_WT_MG-178_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl103_10/0_complex_prep/ENO_AB_WT_MG-010_receptor_2PG_docked.pdb']),

                       ''.join([base_dir_ligand, 'cl223/0_complex_prep/ENO_AB_WT_MG-198_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl223_2/0_complex_prep/ENO_AB_WT_MG-201_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl223_3/0_complex_prep/ENO_AB_WT_MG-210_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl223_4/0_complex_prep/ENO_AB_WT_MG-241_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl223_5/0_complex_prep/ENO_AB_WT_MG-252_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl223_6/0_complex_prep/ENO_AB_WT_MG-195_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl223_7/0_complex_prep/ENO_AB_WT_MG-235_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl223_8/0_complex_prep/ENO_AB_WT_MG-246_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl223_9/0_complex_prep/ENO_AB_WT_MG-250_receptor_2PG_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl223_10/0_complex_prep/ENO_AB_WT_MG-216_receptor_2PG_docked.pdb']),

                       ''.join([base_dir_ligand, 'cl000/0_complex_prep/ENO_AB_WT_MG-001_receptor_2PG_docked.pdb']),
                       ''.join([base_dir, '/ENO_1E9I/0_Crystal_structures/3H8A_A_WT_2PG.pdb'])]

    file_out = ''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/load_2pg_frames.tcl'])
    ligand_list = ['2PG', 'MG']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)

    base_dir_ligand = ''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/1_MG/2_PEP/'])
    file_names_list = [''.join([base_dir_ligand, 'cl302/0_complex_prep/ENO_AB_WT_MG-121_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl302_2/0_complex_prep/ENO_AB_WT_MG-026_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl302_3/0_complex_prep/ENO_AB_WT_MG-065_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl302_4/0_complex_prep/ENO_AB_WT_MG-087_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl302_5/0_complex_prep/ENO_AB_WT_MG-163_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl302_6/0_complex_prep/ENO_AB_WT_MG-109_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl302_7/0_complex_prep/ENO_AB_WT_MG-039_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl302_8/0_complex_prep/ENO_AB_WT_MG-134_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl302_9/0_complex_prep/ENO_AB_WT_MG-154_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl302_10/0_complex_prep/ENO_AB_WT_MG-092_receptor_PEP_docked.pdb']),

                       ''.join([base_dir_ligand, 'cl312/0_complex_prep/ENO_AB_WT_MG-244_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl312_2/0_complex_prep/ENO_AB_WT_MG-099_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl312_3/0_complex_prep/ENO_AB_WT_MG-198_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl312_4/0_complex_prep/ENO_AB_WT_MG-227_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl312_5/0_complex_prep/ENO_AB_WT_MG-251_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl312_6/0_complex_prep/ENO_AB_WT_MG-195_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl312_7/0_complex_prep/ENO_AB_WT_MG-221_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl312_8/0_complex_prep/ENO_AB_WT_MG-231_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl312_9/0_complex_prep/ENO_AB_WT_MG-246_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl312_10/0_complex_prep/ENO_AB_WT_MG-228_receptor_PEP_docked.pdb']),

                       ''.join([base_dir_ligand, 'cl402/0_complex_prep/ENO_AB_WT_MG-058_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl402_2/0_complex_prep/ENO_AB_WT_MG-019_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl402_3/0_complex_prep/ENO_AB_WT_MG-034_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl402_4/0_complex_prep/ENO_AB_WT_MG-057_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl402_5/0_complex_prep/ENO_AB_WT_MG-175_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl402_6/0_complex_prep/ENO_AB_WT_MG-148_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl402_7/0_complex_prep/ENO_AB_WT_MG-113_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl402_8/0_complex_prep/ENO_AB_WT_MG-082_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl402_9/0_complex_prep/ENO_AB_WT_MG-061_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl402_10/0_complex_prep/ENO_AB_WT_MG-041_receptor_PEP_docked.pdb']),

                       ''.join([base_dir_ligand, 'cl313/0_complex_prep/ENO_AB_WT_MG-188_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl313_2/0_complex_prep/ENO_AB_WT_MG-193_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl313_3/0_complex_prep/ENO_AB_WT_MG-209_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl313_4/0_complex_prep/ENO_AB_WT_MG-217_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl313_5/0_complex_prep/ENO_AB_WT_MG-234_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl313_6/0_complex_prep/ENO_AB_WT_MG-098_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl313_7/0_complex_prep/ENO_AB_WT_MG-202_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl313_8/0_complex_prep/ENO_AB_WT_MG-247_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl313_9/0_complex_prep/ENO_AB_WT_MG-176_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl313_10/0_complex_prep/ENO_AB_WT_MG-105_receptor_PEP_docked.pdb']),

                       ''.join([base_dir_ligand, 'cl213/0_complex_prep/ENO_AB_WT_MG-100_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl213_2/0_complex_prep/ENO_AB_WT_MG-101_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl213_3/0_complex_prep/ENO_AB_WT_MG-104_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl213_4/0_complex_prep/ENO_AB_WT_MG-179_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl213_5/0_complex_prep/ENO_AB_WT_MG-230_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl213_6/0_complex_prep/ENO_AB_WT_MG-096_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl213_7/0_complex_prep/ENO_AB_WT_MG-187_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl213_8/0_complex_prep/ENO_AB_WT_MG-192_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl213_9/0_complex_prep/ENO_AB_WT_MG-214_receptor_PEP_docked.pdb']),
                       ''.join([base_dir_ligand, 'cl213_10/0_complex_prep/ENO_AB_WT_MG-181_receptor_PEP_docked.pdb']),

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
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl000_3/0_complex_prep/GAPDH_WT_NAD-131_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl000_4/0_complex_prep/GAPDH_WT_NAD-372_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl000_5/0_complex_prep/GAPDH_WT_NAD-473_receptor_G3P_docked.pdb']),

                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl001/0_complex_prep/GAPDH_WT_NAD-068_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl001_2/0_complex_prep/GAPDH_WT_NAD-343_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl001_3/0_complex_prep/GAPDH_WT_NAD-116_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl001_4/0_complex_prep/GAPDH_WT_NAD-200_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl001_5/0_complex_prep/GAPDH_WT_NAD-550_receptor_G3P_docked.pdb']),

                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl012/0_complex_prep/GAPDH_WT_NAD-159_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl012_2/0_complex_prep/GAPDH_WT_NAD-540_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl012_3/0_complex_prep/GAPDH_WT_NAD-270_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl012_4/0_complex_prep/GAPDH_WT_NAD-344_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl012_5/0_complex_prep/GAPDH_WT_NAD-421_receptor_G3P_docked.pdb']),

                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl010/0_complex_prep/GAPDH_WT_NAD-137_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl010_2/0_complex_prep/GAPDH_WT_NAD-285_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl010_3/0_complex_prep/GAPDH_WT_NAD-380_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl010_4/0_complex_prep/GAPDH_WT_NAD-598_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl010_5/0_complex_prep/GAPDH_WT_NAD-607_receptor_G3P_docked.pdb']),


                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl002/0_complex_prep/GAPDH_WT_NAD-174_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl002_2/0_complex_prep/GAPDH_WT_NAD-214_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl002_3/0_complex_prep/GAPDH_WT_NAD-304_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl002_4/0_complex_prep/GAPDH_WT_NAD-479_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_NAD/1_G3P/cl002_5/0_complex_prep/GAPDH_WT_NAD-542_receptor_G3P_docked.pdb']),
                       ''.join([base_dir, '/GAPDH/0_Crystal_structures/1DC4_chainA.pdb'])]



    file_names_list = [''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl010/0_complex_prep/GAPDH_WT_NAD_rem3PG-275_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl010_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-301_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl010_3/0_complex_prep/GAPDH_WT_NAD_rem3PG-338_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl010_4/0_complex_prep/GAPDH_WT_NAD_rem3PG-358_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl010_5/0_complex_prep/GAPDH_WT_NAD_rem3PG-372_receptor_G3P_docked.pdb']),

                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl100/0_complex_prep/GAPDH_WT_NAD_rem3PG-626_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl100_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-107_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl100_3/0_complex_prep/GAPDH_WT_NAD_rem3PG-007_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl100_4/0_complex_prep/GAPDH_WT_NAD_rem3PG-051_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl100_5/0_complex_prep/GAPDH_WT_NAD_rem3PG-240_receptor_G3P_docked.pdb']),

                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200/0_complex_prep/GAPDH_WT_NAD_rem3PG-025_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-430_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200_3/0_complex_prep/GAPDH_WT_NAD_rem3PG-510_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200_4/0_complex_prep/GAPDH_WT_NAD_rem3PG-612_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200_5/0_complex_prep/GAPDH_WT_NAD_rem3PG-736_receptor_G3P_docked.pdb']),

                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200_6/0_complex_prep/GAPDH_WT_NAD_rem3PG-080_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200_7/0_complex_prep/GAPDH_WT_NAD_rem3PG-140_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200_8/0_complex_prep/GAPDH_WT_NAD_rem3PG-204_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200_9/0_complex_prep/GAPDH_WT_NAD_rem3PG-302_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl200_10/0_complex_prep/GAPDH_WT_NAD_rem3PG-375_receptor_G3P_docked.pdb']),

                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl000/0_complex_prep/GAPDH_WT_NAD_rem3PG-115_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl000_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-306_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl000_3/0_complex_prep/GAPDH_WT_NAD_rem3PG-339_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl000_4/0_complex_prep/GAPDH_WT_NAD_rem3PG-265_receptor_G3P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_NAD_remG3P/1_G3P/cl000_5/0_complex_prep/GAPDH_WT_NAD_rem3PG-004_receptor_G3P_docked.pdb']),

                        ''.join([base_dir, '/GAPDH/0_Crystal_structures/1DC4_chainA.pdb'])]




    file_out = ''.join([base_dir, '/GAPDH/3_MD_post_dock2016/load_g3p_frames.tcl'])
    ligand_list = ['G3P', 'NAD']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)


    base_dir_ligand =''.join([base_dir, '/GAPDH/3_MD_post_dock2016/'])
    file_names_list = [''.join([base_dir_ligand, '1_NAD/2_DPG/cl100/0_complex_prep/GAPDH_WT_NAD-027_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl100_2/0_complex_prep/GAPDH_WT_NAD-139_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl100_3/0_complex_prep/GAPDH_WT_NAD-087_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl100_4/0_complex_prep/GAPDH_WT_NAD-286_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl100_5/0_complex_prep/GAPDH_WT_NAD-359_receptor_13DPG_docked.pdb']),

                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl101/0_complex_prep/GAPDH_WT_NAD-592_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl101_2/0_complex_prep/GAPDH_WT_NAD-234_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl101_3/0_complex_prep/GAPDH_WT_NAD-159_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl101_4/0_complex_prep/GAPDH_WT_NAD-401_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl101_5/0_complex_prep/GAPDH_WT_NAD-491_receptor_13DPG_docked.pdb']),

                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl220/0_complex_prep/GAPDH_WT_NAD-316_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl220_2/0_complex_prep/GAPDH_WT_NAD-399_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl220_3/0_complex_prep/GAPDH_WT_NAD-003_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl220_4/0_complex_prep/GAPDH_WT_NAD-038_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl220_5/0_complex_prep/GAPDH_WT_NAD-146_receptor_13DPG_docked.pdb']),

                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl201/0_complex_prep/GAPDH_WT_NAD-199_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl201_2/0_complex_prep/GAPDH_WT_NAD-307_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl201_3/0_complex_prep/GAPDH_WT_NAD-374_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl201_4/0_complex_prep/GAPDH_WT_NAD-284_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl201_5/0_complex_prep/GAPDH_WT_NAD-514_receptor_13DPG_docked.pdb']),

                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl120/0_complex_prep/GAPDH_WT_NAD-020_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl120_2/0_complex_prep/GAPDH_WT_NAD-071_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl120_3/0_complex_prep/GAPDH_WT_NAD-108_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl120_4/0_complex_prep/GAPDH_WT_NAD-123_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '1_NAD/2_DPG/cl120_5/0_complex_prep/GAPDH_WT_NAD-398_receptor_13DPG_docked.pdb'])]

    file_names_list =[''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl100/0_complex_prep/GAPDH_WT_NAD_rem3PG-010_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl100_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-640_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl100_3/0_complex_prep/GAPDH_WT_NAD_rem3PG-087_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl100_4/0_complex_prep/GAPDH_WT_NAD_rem3PG-236_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl100_5/0_complex_prep/GAPDH_WT_NAD_rem3PG-360_receptor_13DPG_docked.pdb']),

                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl120/0_complex_prep/GAPDH_WT_NAD_rem3PG-242_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl120_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-348_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl120_3/0_complex_prep/GAPDH_WT_NAD_rem3PG-070_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl120_4/0_complex_prep/GAPDH_WT_NAD_rem3PG-141_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl120_5/0_complex_prep/GAPDH_WT_NAD_rem3PG-592_receptor_13DPG_docked.pdb']),

                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl210/0_complex_prep/GAPDH_WT_NAD_rem3PG-516_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl210_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-140_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl210_3/0_complex_prep/GAPDH_WT_NAD_rem3PG-028_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl210_4/0_complex_prep/GAPDH_WT_NAD_rem3PG-422_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl210_5/0_complex_prep/GAPDH_WT_NAD_rem3PG-714_receptor_13DPG_docked.pdb']),

                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl220/0_complex_prep/GAPDH_WT_NAD_rem3PG-035_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl220_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-208_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl220_3/0_complex_prep/GAPDH_WT_NAD_rem3PG-415_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl220_4/0_complex_prep/GAPDH_WT_NAD_rem3PG-591_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl220_5/0_complex_prep/GAPDH_WT_NAD_rem3PG-689_receptor_13DPG_docked.pdb']),

                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl110/0_complex_prep/GAPDH_WT_NAD_rem3PG-176_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl110_2/0_complex_prep/GAPDH_WT_NAD_rem3PG-244_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl110_3/0_complex_prep/GAPDH_WT_NAD_rem3PG-357_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl110_4/0_complex_prep/GAPDH_WT_NAD_rem3PG-641_receptor_13DPG_docked.pdb']),
                       ''.join([base_dir_ligand, '2_NAD_remG3P/2_DPG/cl110_5/0_complex_prep/GAPDH_WT_NAD_rem3PG-737_receptor_13DPG_docked.pdb']),

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
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl101_3/0_complex_prep/TALB_WT_APO-244_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl101_4/0_complex_prep/TALB_WT_APO-388_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl101_5/0_complex_prep/TALB_WT_APO-470_receptor_S7P_docked.pdb']),

                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl102/0_complex_prep/TALB_WT_APO-091_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl102_2/0_complex_prep/TALB_WT_APO-535_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl102_3/0_complex_prep/TALB_WT_APO-282_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl102_4/0_complex_prep/TALB_WT_APO-433_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl102_5/0_complex_prep/TALB_WT_APO-752_receptor_S7P_docked.pdb']),

                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl212/0_complex_prep/TALB_WT_APO-069_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl212_2/0_complex_prep/TALB_WT_APO-777_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl212_3/0_complex_prep/TALB_WT_APO-190_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl212_4/0_complex_prep/TALB_WT_APO-357_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl212_5/0_complex_prep/TALB_WT_APO-650_receptor_S7P_docked.pdb']),

                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl213/0_complex_prep/TALB_WT_APO-030_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl213_2/0_complex_prep/TALB_WT_APO-732_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl213_3/0_complex_prep/TALB_WT_APO-675_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl213_4/0_complex_prep/TALB_WT_APO-720_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl213_5/0_complex_prep/TALB_WT_APO-755_receptor_S7P_docked.pdb']),

                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl211/0_complex_prep/TALB_WT_APO-089_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl211_2/0_complex_prep/TALB_WT_APO-155_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl211_3/0_complex_prep/TALB_WT_APO-236_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl211_4/0_complex_prep/TALB_WT_APO-383_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '1_APO/1_S7P_linear/cl211_5/0_complex_prep/TALB_WT_APO-400_receptor_S7P_docked.pdb'])]

    file_names_list = [''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl101/0_complex_prep/TALB_WT_HALO_S7P_remS7P-237_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl101_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-260_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl101_3/0_complex_prep/TALB_WT_HALO_S7P_remS7P-305_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl101_4/0_complex_prep/TALB_WT_HALO_S7P_remS7P-177_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl101_5/0_complex_prep/TALB_WT_HALO_S7P_remS7P-333_receptor_S7P_linear_docked.pdb']),

                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl302/0_complex_prep/TALB_WT_HALO_S7P_remS7P-091_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl302_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-492_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl302_3/0_complex_prep/TALB_WT_HALO_S7P_remS7P-235_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl302_4/0_complex_prep/TALB_WT_HALO_S7P_remS7P-287_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl302_5/0_complex_prep/TALB_WT_HALO_S7P_remS7P-309_receptor_S7P_linear_docked.pdb']),

                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl322/0_complex_prep/TALB_WT_HALO_S7P_remS7P-042_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl322_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-184_receptor_S7P_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl322_3/0_complex_prep/TALB_WT_HALO_S7P_remS7P-148_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl322_4/0_complex_prep/TALB_WT_HALO_S7P_remS7P-225_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl322_5/0_complex_prep/TALB_WT_HALO_S7P_remS7P-273_receptor_S7P_linear_docked.pdb']),

                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl131/0_complex_prep/TALB_WT_HALO_S7P_remS7P-114_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl131_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-140_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl131_3/0_complex_prep/TALB_WT_HALO_S7P_remS7P-207_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl131_4/0_complex_prep/TALB_WT_HALO_S7P_remS7P-051_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl131_5/0_complex_prep/TALB_WT_HALO_S7P_remS7P-400_receptor_S7P_linear_docked.pdb']),

                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl122/0_complex_prep/TALB_WT_HALO_S7P_remS7P-028_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl122_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-029_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl122_3/0_complex_prep/TALB_WT_HALO_S7P_remS7P-351_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl122_4/0_complex_prep/TALB_WT_HALO_S7P_remS7P-187_receptor_S7P_linear_docked.pdb']),
                        ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/1_S7P_linear/cl122_5/0_complex_prep/TALB_WT_HALO_S7P_remS7P-424_receptor_S7P_linear_docked.pdb']),

                       ''.join([base_dir, '/TalB/0_Crystalline_structures/1UCW_wt_S7P_watbox.pdb'])]

    file_out = ''.join([base_dir, '/TalB/3_MD_post_dock2016/load_s7p_halo_frames.tcl'])
    ligand_list = ['S7P']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)


    base_dir_ligand = ''.join([base_dir, '/TalB/3_MD_post_dock2016/'])
    file_names_list = [''.join([base_dir_ligand, '1_APO/2_F6P/cl202/0_complex_prep/TALB_WT_APO-533_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl202_2/0_complex_prep/TALB_WT_APO-314_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl202_3/0_complex_prep/TALB_WT_APO-176_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl202_4/0_complex_prep/TALB_WT_APO-440_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl202_5/0_complex_prep/TALB_WT_APO-297_receptor_F6P_docked.pdb']),

                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl212/0_complex_prep/TALB_WT_APO-139_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl212_2/0_complex_prep/TALB_WT_APO-464_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl212_3/0_complex_prep/TALB_WT_APO-091_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl212_4/0_complex_prep/TALB_WT_APO-287_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl212_5/0_complex_prep/TALB_WT_APO-193_receptor_F6P_docked.pdb']),

                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl312/0_complex_prep/TALB_WT_APO-799_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl312_2/0_complex_prep/TALB_WT_APO-003_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl312_3/0_complex_prep/TALB_WT_APO-203_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl312_4/0_complex_prep/TALB_WT_APO-504_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl312_5/0_complex_prep/TALB_WT_APO-669_receptor_F6P_docked.pdb']),

                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl201/0_complex_prep/TALB_WT_APO-179_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl201_2/0_complex_prep/TALB_WT_APO-283_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl201_3/0_complex_prep/TALB_WT_APO-341_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl201_4/0_complex_prep/TALB_WT_APO-396_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl201_5/0_complex_prep/TALB_WT_APO-435_receptor_F6P_docked.pdb']),

                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl101/0_complex_prep/TALB_WT_APO-373_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl101_2/0_complex_prep/TALB_WT_APO-385_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl101_3/0_complex_prep/TALB_WT_APO-404_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl101_4/0_complex_prep/TALB_WT_APO-166_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '1_APO/2_F6P/cl101_5/0_complex_prep/TALB_WT_APO-368_receptor_F6P_docked.pdb'])]

    file_names_list = [''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl101/0_complex_prep/TALB_WT_HALO_S7P_remS7P-492_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl101_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-056_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl101_3/0_complex_prep/TALB_WT_HALO_S7P_remS7P-159_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl101_4/0_complex_prep/TALB_WT_HALO_S7P_remS7P-245_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl101_5/0_complex_prep/TALB_WT_HALO_S7P_remS7P-345_receptor_F6P_docked.pdb']),

                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl201/0_complex_prep/TALB_WT_HALO_S7P_remS7P-381_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl201_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-117_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl201_3/0_complex_prep/TALB_WT_HALO_S7P_remS7P-051_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl201_4/0_complex_prep/TALB_WT_HALO_S7P_remS7P-200_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl201_5/0_complex_prep/TALB_WT_HALO_S7P_remS7P-292_receptor_F6P_docked.pdb']),


                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl311/0_complex_prep/TALB_WT_HALO_S7P_remS7P-014_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl311_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-222_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl311_3/0_complex_prep/TALB_WT_HALO_S7P_remS7P-090_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl311_4/0_complex_prep/TALB_WT_HALO_S7P_remS7P-181_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl311_5/0_complex_prep/TALB_WT_HALO_S7P_remS7P-376_receptor_F6P_docked.pdb']),

                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl301/0_complex_prep/TALB_WT_HALO_S7P_remS7P-028_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl301_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-102_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl301_3/0_complex_prep/TALB_WT_HALO_S7P_remS7P-193_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl301_4/0_complex_prep/TALB_WT_HALO_S7P_remS7P-336_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl301_5/0_complex_prep/TALB_WT_HALO_S7P_remS7P-414_receptor_F6P_docked.pdb']),

                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl211/0_complex_prep/TALB_WT_HALO_S7P_remS7P-074_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl211_2/0_complex_prep/TALB_WT_HALO_S7P_remS7P-107_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl211_3/0_complex_prep/TALB_WT_HALO_S7P_remS7P-216_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl211_4/0_complex_prep/TALB_WT_HALO_S7P_remS7P-274_receptor_F6P_docked.pdb']),
                       ''.join([base_dir_ligand, '2_HALO_S7P_remS7P/2_F6P/cl211_5/0_complex_prep/TALB_WT_HALO_S7P_remS7P-304_receptor_F6P_docked.pdb']),

                       ''.join([base_dir, '/TalB/0_Crystalline_structures/1UCW_wt_S7P_watbox.pdb'])]


    file_out = ''.join([base_dir, '/TalB/3_MD_post_dock2016/load_f6p_halo_frames.tcl'])
    ligand_list = ['F6P']
    create_vmd_file(file_names_list, file_out, ligand_list, binding_residues)


def main():
    base_dir = '/home/mrama/Desktop/MD/eMASS-MD_complete_data/MD_data/'
    #create_tcl_for_eno_1e9i(base_dir)
    #create_tcl_for_eno_ab(base_dir)
    #create_tcl_for_gapd(base_dir)
    create_tcl_for_talb(base_dir)

if __name__ == '__main__':
    main()

