import os
import subprocess
from shutil import copy, copyfile, rmtree

import jinja2
import yaml


def discard_first_n_frames_from_traj(output_dir, md_traj_dir, prmtop_file, md_traj_in, md_traj_out, n_frames):
    """
        Takes a MD trajectory, removes the first n frames and writes the resulting trajectory.

    :param output_dir: path to the folder where the cpptraj files will be written
    :param md_traj_dir: path to where the MD simulation is
    :param prmtop_file: prmtop file path
    :param md_traj_in: input trajectory name
    :param md_traj_out: output trajectory name
    :param n_frames: integer specifying how many initial frames to discard
    :return: None
    """

    with open(''.join([output_dir, '/discard_n_frames.in']), 'w') as f_cpp_complex_solv:
        f_cpp_complex_solv.write(' '.join(['trajin', md_traj_in, str(n_frames), 'last 1\n']) )
        f_cpp_complex_solv.write(''.join(['trajout ', md_traj_out, ' crd\n']))
        f_cpp_complex_solv.write('go')

    with open('/'.join([output_dir, 'run_cpptraj.sh']), 'w') as f_run_cpptraj:
        f_run_cpptraj.write(''.join(['cd ', md_traj_dir, '\n']))
        f_run_cpptraj.write(
            ''.join(['$AMBERHOME/bin/cpptraj -i ', output_dir, '/discard_n_frames.in -p ', prmtop_file, '\n']))

    os.chmod('./run_cpptraj.sh', 0755)
    subprocess.call('./run_cpptraj.sh', shell=True)


# TODO: make it based on pytraj or something python
def prep_cpptraj_files(output_dir, md_traj_name, ligand, cpp_names, frame=None):
    """
        Creates 4 cpptraj scripts to generate the following pdb files from a given MD trajectory:
        complex_solv.pdb, complex_vac.pdb, ligand_vac, and receptor_vac.

    :param output_dir: path to the folder where the cpptraj files will be written
    :param md_traj_name: base name of the MD simulation, e.g. G6PD_WT_NADP_6PGL_docked
    :param cl: the cluster name, e.g. cl101
    :param ligand: name or number of the ligand bound to the enzyme, as in the pdb
    :param cpp_names: a list with the names of the pdb files to be created
    :param md_sim: if is not None it specifies the number of the MD simulation.
    :param frame: specify from which trajectory frame will the pdb files be created, if None it gets the last frame
    :return: None
    """

    # write cpptraj script for solvated complex
    with open(''.join([output_dir, '/get_traj_frame_', cpp_names[0], '.in']), 'w') as f_cpp_complex_solv:
        trajin_com = ' '.join(['trajin', md_traj_name, str(frame), str(frame), '\n']) if frame else ''.join(['trajin ', md_traj_name, ' lastframe\n'])
        f_cpp_complex_solv.write(trajin_com)
        f_cpp_complex_solv.write(''.join(['trajout ', cpp_names[0], '.pdb pdb\n']))
        f_cpp_complex_solv.write('go')

    # write cpptraj script for the complex in vacuum
    with open(''.join([output_dir, '/get_traj_frame_', cpp_names[1], '.in']), 'w') as f_cpp_complex_vac:
        trajin_com = ' '.join(['trajin', md_traj_name, str(frame), str(frame), '\n']) if frame else ''.join(['trajin ', md_traj_name, ' lastframe\n'])
        f_cpp_complex_vac.write(trajin_com)
        f_cpp_complex_vac.write(''.join(['trajout ', cpp_names[1], '.pdb pdb\n']))
        f_cpp_complex_vac.write('strip :WAT\n')
        f_cpp_complex_vac.write('strip :Na+\n')
        f_cpp_complex_vac.write('go')

    # write cpptraj script for the ligand in vacuum
    with open(''.join([output_dir, '/get_traj_frame_', cpp_names[2], '.in']), 'w') as f_cpp_ligand_vac:
        trajin_com = ' '.join(['trajin', md_traj_name, str(frame), str(frame), '\n']) if frame else ''.join(['trajin ', md_traj_name, ' lastframe\n'])
        f_cpp_ligand_vac.write(trajin_com)
        f_cpp_ligand_vac.write(''.join(['trajout ', cpp_names[2], '.pdb pdb\n']))
        f_cpp_ligand_vac.write(''.join(['strip !(:', ligand, ')\n']))
        f_cpp_ligand_vac.write('go')

    # write cpptraj script for receptor in vacuum
    with open(''.join([output_dir, '/get_traj_frame_', cpp_names[3], '.in']), 'w') as f_cpp_receptor_vac:
        trajin_com = ' '.join(['trajin', md_traj_name, str(frame), str(frame), '\n']) if frame else ''.join(['trajin ', md_traj_name, ' lastframe\n'])
        f_cpp_receptor_vac.write(trajin_com)
        f_cpp_receptor_vac.write(''.join(['trajout ', cpp_names[3], '.pdb pdb\n']))
        f_cpp_receptor_vac.write('strip :WAT\n')
        f_cpp_receptor_vac.write('strip :Na+\n')
        f_cpp_receptor_vac.write(''.join(['strip :', ligand, '\n']))
        f_cpp_receptor_vac.write('go')


def prep_get_frame_from_traj(output_dir, md_traj_dir, prmtop_file, cpp_names):
    """
        Creates the bash script that will run the cpptraj files that actually create pdb files from a given MD
        trajectory for the complex of enzyme+ligand in water and vacuum, the receptor in vacuum, and the ligand in
        vacuum.

    :param output_dir: path to where the bash script to run cpptraj will be written, and where the cpptraj scripts are
    :param md_traj_dir: path to where the MD simulation is
    :param prmtop_file: prmtop file path
    :param cpp_files: a list with the names of the pdb files to be created
    :return: None
    """

    # write bash script to run cpptraj
    with open('/'.join([output_dir, 'run_cpptraj_mmpbsa.sh']), 'w') as f_run_cpptraj:
        f_run_cpptraj.write(''.join(['cd ', md_traj_dir, '\n']))
        f_run_cpptraj.write(
            ''.join(['$AMBERHOME/bin/cpptraj -i ', output_dir, '/get_traj_frame_', cpp_names[0], '.in -p ', prmtop_file, '\n']))
        f_run_cpptraj.write(
            ''.join(['$AMBERHOME/bin/cpptraj -i ', output_dir, '/get_traj_frame_', cpp_names[1], '.in -p ', prmtop_file, '\n']))
        f_run_cpptraj.write(
            ''.join(['$AMBERHOME/bin/cpptraj -i ', output_dir, '/get_traj_frame_', cpp_names[2], '.in -p ', prmtop_file, '\n']))
        f_run_cpptraj.write(
            ''.join(['$AMBERHOME/bin/cpptraj -i ', output_dir, '/get_traj_frame_', cpp_names[3], '.in -p ', prmtop_file, '\n']))


# TODO: write unit test for case with multiple ligands (take ENO_AB with MG, 2PG)
def prep_complex_pdbs(mmpbsa_base_dir, mmpbsa_dir, md_traj_dir, cl, pdb_names, leaprc_base, ligand_params, ligand,
                      frame=None):

    """
        This function copies all the pdb files for complex_solv, complex_vac, receptor_vac, and ligand_vac to the
        0_complex_prep folder of the respective cluster (cl), and creates the leaprc files that are used to
        create the respective inpcrd and prmtop files. In the end, tleap is executed.

    :param mmpbsa_base_dir: directory where MMPBSA is, excluding cluster info, e.g. /home/marta/G6PD/4_MMPBSA/1_APO/1_G6P
    :param mmpbsa_dir: directory where MMPBSA is
    :param md_traj_dir: absolute path to where the MD simulation is, and thus where the PDB files are
    :param cl: the cluster name, e.g. cl101
    :param pdb_names: a list with the names of the pdb files to be created
    :param leaprc_base: absolute path to where the base for the leaprc files is located
    :param ligand_params: absolute path to where the ligand parameters are, or a list of these
    :param ligand: the name of the ligand bound to the enzyme, or a list of these
    :param frame: specify frame to be extracted if doing mmpbsa per frame
    :return: None
    """

    # create folders for MMPBSA
    complex_prep_dir = ''.join([mmpbsa_dir, '/frame_', str(frame), '/0_complex_prep']) if frame else '/'.join([mmpbsa_base_dir, cl, '0_complex_prep'])
    if not os.path.exists(complex_prep_dir):
        os.makedirs(complex_prep_dir)
    os.chdir(complex_prep_dir)

    for pdb_name in pdb_names:
        # copy pdb files
        copy('/'.join([md_traj_dir, pdb_name + '.pdb']), './')

        # create leaprc file
        copyfile(leaprc_base, ''.join(['./leaprc_', pdb_name]))

        with open(''.join(['leaprc_' + pdb_name]), 'a') as f_leaprc:
            if type(ligand_params) is list:
                for i, param_set in enumerate(ligand_params):
                    f_leaprc.write(''.join(['loadamberparams ', param_set, '/', ligand[i], '.frcmod\n']))
                    f_leaprc.write(''.join(['loadOff ', param_set, '/', ligand[i], '.off\n']))
            else:
                f_leaprc.write(''.join(['loadamberparams ', ligand_params, '/', ligand, '.frcmod\n']))
                f_leaprc.write(''.join(['loadOff ', ligand_params, '/', ligand, '.off\n']))

            f_leaprc.write('set default PBRadii bondi\n')
            f_leaprc.write(''.join(['mol = loadpdb ',  pdb_name, '.pdb\n']))
            f_leaprc.write(''.join(['saveamberparm mol ', pdb_name, '.prmtop ', pdb_name, '.inpcrd\n']))
            f_leaprc.write('quit')

        # run leaprc file
        subprocess.call(''.join(['tleap -f leaprc_', pdb_name]), shell=True)


def prep_mmpbsa_run(mmpbsa_dir, mmpbsa_in_file, md_traj, frame=None):
    """
        Creates the script to run MMPBSA, after copying the mmpbsa_pb.in file to the same folder.

    :param mmpbsa_dir: path to where the script is located and from which MMPBSA will be executed.
    :param mmpbsa_in_file: path to the mmpbsa_pb.in file
    :param md_traj: path to the MD simulation used in the MMPBSA
    :param frame:
    :return: None
    """

    if not frame:
        copy(mmpbsa_in_file, ''.join([mmpbsa_dir, '/']))

    complex_solv = '0_complex_prep/complex_solv.prmtop'
    complex_vac = '0_complex_prep/complex_vac.prmtop'
    receptor_vac = '0_complex_prep/receptor_vac.prmtop'
    ligand_vac = '0_complex_prep/ligand_vac.prmtop'

    frame_dir = ''.join(['frame_', str(frame)])
    mmpbsa_dir = mmpbsa_dir if not frame else '/'.join([mmpbsa_dir, frame_dir])

    with open('/'.join([mmpbsa_dir, 'run_mmpbsa.sh']), 'w') as f_run_mmpbsa:
        f_run_mmpbsa.write(''.join(['$AMBERHOME/bin/MMPBSA.py -O -i mmpbsa_pb.in -o FINAL_RESULTS_MMPBSA.dat -sp ',
                                    complex_solv, ' -cp ', complex_vac, ' -rp ', receptor_vac, ' -lp ', ligand_vac,
                                    ' -y ', md_traj]))

    os.chmod('/'.join([mmpbsa_dir, 'run_mmpbsa.sh']), 0755)


def prep_mmpbsa_in(template_folder, template_file, mmpbsa_dict):
    """
        Creates MMPBSA input script (to specify gb, pb, number of frames, other options)
        Requires template file: 'mmpbsa_perframe.jinja'
        Returns filename of input script

    :param template_folder:
    :param template_file:
    :param mmpbsa_dict:
    :return: None
    """

    # allows loading of templates off filesystems
    templateEnv = jinja2.Environment(loader=jinja2.FileSystemLoader(template_folder), trim_blocks=True)

    # load the template
    template = templateEnv.get_template(template_file)

    # process the template to produce our final text.
    outputText = template.render(mmpbsa_dict)

    file_out = mmpbsa_dict['filename']
    with open(file_out, "w") as f:
        f.write(outputText)




# options to use for running MMPBSA perframe (leave this the same)
pb_perframe = yaml.load("""
---
gen_options:
    - keep_files=0
    - verbose=2
pb_options:
    - fillratio=4.0
    - inp=1
    - radiopt=0
""")


def mmpbsa_per_frame_eno():

    per_frame = True
    reorder_atoms = False
    discard_n_frames = 7

    enzyme = 'ENO_1E9I'
    form = 'APO'
    form_folder = '_'.join(['1', form])
    ligand = '2PG'
    ligand_folder = '_'.join(['1', ligand])

    #cl_frames = OrderedDict((('cl101', 70), ('cl302', 70), ('cl322', 140)))  # 725 G6PD
    cl_frames = {'cl200_4': 210}  # 725 G6PD
    # md_sim_list = ['1']
    md_sim_list = ['1']

    base_form = '_'.join(['1E9I_WT_APO', ligand, 'docked'])
    base_dir = '/'.join(['/home/marta', enzyme])
    scripts_dir = '/home/marta/scripts/mmpbsa_scripts'
    cpp_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']
    traj_format_list = ['crd']
    traj_format_i = 0
    traj_type = '100'

    for i, (cl, n_frames) in enumerate(cl_frames.items()):
        traj_format = traj_format_list[i]
        md_traj_name = ''.join([base_form, '_', cl, '_', md_sim_list[i], '_', traj_type, '.', traj_format]) if md_sim_list[i] else \
                                ''.join([base_form, '_', cl, '_', traj_type, '.', traj_format])

        md_traj_dir = '/'.join([base_dir, '3_MD_post_dock2016', form_folder, ligand_folder, cl])
        md_traj = '/'.join([md_traj_dir, md_traj_name])

        prmtop_file = '.'.join([base_form, 'prmtop'])
        leaprc_base = '/'.join([scripts_dir, 'base_files/base_leaprc'])
        mmpbsa_in_file = '/'.join([scripts_dir, 'base_files/mmpbsa_pb.in'])

        mmpbsa_base_dir = '/'.join([base_dir, '4_MMPBSA2016', form_folder, ligand_folder])

        if per_frame:
            mmpbsa_dir = '/'.join([mmpbsa_base_dir,  cl, 'per_frame'])
            if not os.path.exists(mmpbsa_dir):
                os.makedirs(mmpbsa_dir)
            prep_get_frame_from_traj(mmpbsa_dir, md_traj_dir, prmtop_file, cpp_names)

            for frame in range(1, n_frames+1):
                os.chdir(mmpbsa_dir)
                frame_dir = ''.join([mmpbsa_dir, '/frame_', str(frame)])

                # generate cpptraj files
                #ligand = 'DPG,NAD'
                prep_cpptraj_files(mmpbsa_dir, md_traj_name, ligand, cpp_names, frame)

                # generate pdb files for complex_solv/vac receptor/ligand_vac
                os.chdir(mmpbsa_dir)
                os.chmod('./run_cpptraj_mmpbsa.sh', 0755)
                subprocess.call('./run_cpptraj_mmpbsa.sh', shell=True)

                # prep the complex for mmpbsa
                ligands = [ligand, 'MG']
                #ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '2_PEP_liz']),
                ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '1_2PG']),
                                      '/'.join([base_dir, '0_Ligands_param', '3_MG'])]
                #ligand_params = '/'.join([base_dir, '0_Ligands_param', '3_S7P_linear_noopt'])

                prep_complex_pdbs(mmpbsa_base_dir, mmpbsa_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params_list, ligands, frame)
                #prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params, ligand, frame)

                # generate an mmpbsa.in file for the current frame
                mmpbsa_dict = pb_perframe
                mmpbsa_dict['filename'] = frame_dir + '/mmpbsa_pb.in'
                mmpbsa_dict['start_frame'] = frame
                mmpbsa_dict['end_frame'] = frame

                template_folder = '/'.join([scripts_dir, 'base_files'])
                mmpbsa_per_frame_template = 'mmpbsa_perframe.jinja'
                prep_mmpbsa_in(template_folder, mmpbsa_per_frame_template, mmpbsa_dict)

                # prepare script to run mmpbsa for the current frame
                prep_mmpbsa_run(mmpbsa_dir, mmpbsa_in_file, md_traj, frame)

                # run mmpbsa for the current frame
                os.chdir(frame_dir)
                os.chmod('run_mmpbsa.sh', 0755)
                subprocess.call('./run_mmpbsa.sh', shell=True)

                # remove all files except for the results file
                try:
                    os.remove('_*')
                except OSError:
                    pass

                os.remove('mmpbsa_pb.in')
                os.remove('run_mmpbsa.sh')
                rmtree('0_complex_prep')
                # mmpbsa takes about 1min
                #time.sleep(15)

        else:

            mmpbsa_dir = '/'.join([mmpbsa_base_dir,  cl])
            os.chdir(scripts_dir)

            if discard_n_frames:
                print md_traj
                print md_traj_dir
                md_traj_name_out = ''.join([md_traj_name[:-4], '_discard7f.crd'])
                discard_first_n_frames_from_traj(scripts_dir, md_traj_dir, prmtop_file, md_traj_name, md_traj_name_out, discard_n_frames)
                md_traj_name = md_traj_name_out
                md_traj = '/'.join([md_traj_dir, md_traj_name])
                print md_traj
                print '---'


            ligand = 'PEP'
            prep_cpptraj_files(scripts_dir, md_traj_name, ligand, cpp_names)
            prep_get_frame_from_traj(scripts_dir, md_traj_dir, prmtop_file, cpp_names)  # remove?


            os.chmod('./run_cpptraj_mmpbsa.sh', 0755)
            subprocess.call('./run_cpptraj_mmpbsa.sh', shell=True)

            # if reorder_atoms:
            #
            #     pdb_base = '/'.join([base_dir, '4_MMPBSA', form_folder, ligand_folder, cl, '0_complex_prep'])
            #     off_file = '/'.join([base_dir, '0_Ligands_param', ligand_folder, ligand])
            #     n_atoms = 15
            #
            #     for cpp_name in cpp_names:
            #         pdb_file = '/'.join([pdb_base, cpp_name])
            #         reorder_ligand_atoms(pdb_file, off_file, ligand, n_atoms)


            #ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '1_G3P'])],
            #                          '/'.join([base_dir, '0_Ligands_param', '3_NAD']),
            #                          '/'.join([base_dir, '0_Ligands_param', '8_H2PO4_opt'])]

            ligands = ['MG', 'PEP']
            ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '3_MG']),
                                  '/'.join([base_dir, '0_Ligands_param', '2_PEP_liz'])]
            prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params_list, ligands)

            #ligand_params = '/'.join([base_dir, '0_Ligands_param', '2_F6P'])
            #prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params, ligand)
            prep_mmpbsa_run(mmpbsa_dir, mmpbsa_in_file, md_traj)
            os.chdir(mmpbsa_dir)

            # run mmpbsa for the current frame
            os.chmod('./run_mmpbsa.sh', 0755)
            subprocess.call('nohup ./run_mmpbsa.sh &', shell=True)

        traj_format_i += 1


def mmpbsa_per_frame_eno_ab():

    per_frame = True
    reorder_atoms = False
    discard_n_frames = 7

    enzyme = 'ENO_AB'
    form = 'MG'
    form_folder = '_'.join(['1', form])
    ligand = 'PEP'
    ligand_folder = '_'.join(['2', ligand])

    #cl_frames = OrderedDict((('cl101', 70), ('cl302', 70), ('cl322', 140)))  # 725 G6PD
    cl_frames = {#'cl302_11': 70, 'cl302_12': 70, 'cl302_13': 70, 'cl302_14': 70, 'cl302_15': 70}
                #'cl302_16': 70, 'cl302_17': 70, 'cl302_19': 70, 'cl302_20': 70,
                #'cl202': 70} #, 'cl302_18': 70,
                #'cl202_2': 70, 'cl202_3': 70, 'cl202_4': 70, 'cl202_5': 70,
                #'cl202_6': 70} #, 
                'cl302_16': 70, 'cl302_17': 70, 'cl302_19': 70, 'cl302_20': 70}

    """cl_frames = {'cl302': 70, 'cl302_2': 70, 'cl302_3': 70, 'cl302_4': 70, 'cl302_5': 70,
                 'cl302_6': 70, 'cl302_7': 70, 'cl302_8': 70, 'cl302_9': 70,
                 'cl312': 70, 'cl312_2': 70, 'cl312_3': 70, 'cl312_4': 70, 'cl312_5': 70,
                 'cl312_6': 70, 'cl312_7': 70, 'cl312_8': 70, 'cl312_9': 70, 'cl312_10': 70,
                 'cl402': 70, 'cl402_2': 70, 'cl402_3': 70, 'cl402_4': 70, 'cl402_5': 70,
                 'cl402_6': 70, 'cl402_7': 70, 'cl402_8': 70, 'cl402_9': 70, 'cl402_10': 70,
                 'cl313': 70, 'cl313_2': 70, 'cl313_3': 70, 'cl313_4': 70, 'cl313_5': 70,
                 'cl313_6': 70, 'cl313_7': 70, 'cl313_8': 70, 'cl313_9': 70, 'cl313_10': 70}"""

    # md_sim_list = ['1']
    md_sim_list = ['1' for i in range(len(cl_frames))]

    base_form = '_'.join(['ENO_AB_WT_APO', ligand, 'docked'])
    base_dir = '/'.join(['/home/marta', enzyme])
    scripts_dir = '/home/marta/scripts/mmpbsa_scripts'
    cpp_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']
    traj_format_list = ['crd' for i in range(len(cl_frames))]
    traj_format_i = 0
    traj_type = '10'

    for i, (cl, n_frames) in enumerate(cl_frames.items()):
        traj_format = traj_format_list[i]
        md_traj_name = ''.join([base_form, '_', cl, '_', md_sim_list[i], '_', traj_type, '.', traj_format]) if md_sim_list[i] else \
                                ''.join([base_form, '_', cl, '_', traj_type, '.', traj_format])

        md_traj_dir = '/'.join([base_dir, '3_MD_post_dock2016', form_folder, ligand_folder, cl])
        md_traj = '/'.join([md_traj_dir, md_traj_name])

        prmtop_file = '.'.join([base_form, 'prmtop'])
        leaprc_base = '/'.join([scripts_dir, 'base_files/base_leaprc'])
        mmpbsa_in_file = '/'.join([scripts_dir, 'base_files/mmpbsa_pb.in'])

        mmpbsa_base_dir = '/'.join([base_dir, '4_MMPBSA2016', form_folder, ligand_folder])

        if per_frame:
            mmpbsa_dir = '/'.join([mmpbsa_base_dir,  cl, 'per_frame'])
            if not os.path.exists(mmpbsa_dir):
                os.makedirs(mmpbsa_dir)
            prep_get_frame_from_traj(mmpbsa_dir, md_traj_dir, prmtop_file, cpp_names)

            for frame in range(1, n_frames+1):
                os.chdir(mmpbsa_dir)
                frame_dir = ''.join([mmpbsa_dir, '/frame_', str(frame)])

                # generate cpptraj files
                prep_cpptraj_files(mmpbsa_dir, md_traj_name, ligand, cpp_names, frame)

                # generate pdb files for complex_solv/vac receptor/ligand_vac
                os.chdir(mmpbsa_dir)
                os.chmod('./run_cpptraj_mmpbsa.sh', 0755)
                subprocess.call('./run_cpptraj_mmpbsa.sh', shell=True)

                # prep the complex for mmpbsa
                ligands = [ligand, 'MG']
                ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '2_PEP_liz']),
                #ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '7_2PG_liz']),
                                      '/'.join([base_dir, '0_Ligands_param', '3_MG'])]
                #ligand_params = '/'.join([base_dir, '0_Ligands_param', '3_S7P_linear_noopt'])

                prep_complex_pdbs(mmpbsa_base_dir, mmpbsa_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params_list, ligands, frame)
                #prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params, ligand, frame)

                # generate an mmpbsa.in file for the current frame
                mmpbsa_dict = pb_perframe
                mmpbsa_dict['filename'] = frame_dir + '/mmpbsa_pb.in'
                mmpbsa_dict['start_frame'] = frame
                mmpbsa_dict['end_frame'] = frame

                template_folder = '/'.join([scripts_dir, 'base_files'])
                mmpbsa_per_frame_template = 'mmpbsa_perframe.jinja'
                prep_mmpbsa_in(template_folder, mmpbsa_per_frame_template, mmpbsa_dict)

                # prepare script to run mmpbsa for the current frame
                prep_mmpbsa_run(mmpbsa_dir, mmpbsa_in_file, md_traj, frame)

                # run mmpbsa for the current frame
                os.chdir(frame_dir)
                os.chmod('run_mmpbsa.sh', 0755)
                subprocess.call('./run_mmpbsa.sh', shell=True)

                # remove all files except for the results file
                try:
                    os.remove('_*')
                except OSError:
                    pass

                os.remove('mmpbsa_pb.in')
                os.remove('run_mmpbsa.sh')
                rmtree('0_complex_prep')
                # mmpbsa takes about 1min
                #time.sleep(15)
        """
        else:

            mmpbsa_dir = '/'.join([mmpbsa_base_dir,  cl])
            os.chdir(scripts_dir)

            if discard_n_frames:
                print md_traj
                print md_traj_dir
                md_traj_name_out = ''.join([md_traj_name[:-4], '_discard7f.crd'])
                discard_first_n_frames_from_traj(scripts_dir, md_traj_dir, prmtop_file, md_traj_name, md_traj_name_out, discard_n_frames)
                md_traj_name = md_traj_name_out
                md_traj = '/'.join([md_traj_dir, md_traj_name])
                print md_traj
                print '---'


            ligand = 'PEP'
            prep_cpptraj_files(scripts_dir, md_traj_name, ligand, cpp_names)
            prep_get_frame_from_traj(scripts_dir, md_traj_dir, prmtop_file, cpp_names)  # remove?


            os.chmod('./run_cpptraj_mmpbsa.sh', 0755)
            subprocess.call('./run_cpptraj_mmpbsa.sh', shell=True)

            # if reorder_atoms:
            #
            #     pdb_base = '/'.join([base_dir, '4_MMPBSA', form_folder, ligand_folder, cl, '0_complex_prep'])
            #     off_file = '/'.join([base_dir, '0_Ligands_param', ligand_folder, ligand])
            #     n_atoms = 15
            #
            #     for cpp_name in cpp_names:
            #         pdb_file = '/'.join([pdb_base, cpp_name])
            #         reorder_ligand_atoms(pdb_file, off_file, ligand, n_atoms)


            #ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '1_G3P'])],
            #                          '/'.join([base_dir, '0_Ligands_param', '3_NAD']),
            #                          '/'.join([base_dir, '0_Ligands_param', '8_H2PO4_opt'])]

            ligands = ['MG', 'PEP']
            ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '3_MG']),
                                  '/'.join([base_dir, '0_Ligands_param', '2_PEP_liz'])]
            prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params_list, ligands)

            #ligand_params = '/'.join([base_dir, '0_Ligands_param', '2_F6P'])
            #prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params, ligand)
            prep_mmpbsa_run(mmpbsa_dir, mmpbsa_in_file, md_traj)
            os.chdir(mmpbsa_dir)

            # run mmpbsa for the current frame
            os.chmod('./run_mmpbsa.sh', 0755)
            subprocess.call('nohup ./run_mmpbsa.sh &', shell=True)
        """
        traj_format_i += 1


def mmpbsa_per_frame_gapdh():

    per_frame = True
    reorder_atoms = False
    discard_n_frames = 7

    enzyme = 'GAPDH'
    form = 'NAD'
    form_folder = '_'.join(['2', form])
    ligand = '13DPG'
    ligand_folder = '_'.join(['2', ligand])

    cl_frames = {'cl201': 70, 'cl201_2': 70, 'cl201_3': 70, 'cl201_4': 70, 'cl201_5': 70}
    #cl_frames = {'cl001': 70, 'cl001_2': 70, 'cl001_3': 70, 'cl001_4': 70, 'cl001_5': 70,
    #             'cl012': 70, 'cl012_2': 70, 'cl012_3': 70, 'cl012_4': 70, 'cl012_5': 70,
    #             'cl010': 70, 'cl010_2': 70, 'cl010_3': 70, 'cl010_4': 70, 'cl010_5': 70,
    #             'cl000': 70, 'cl000_2': 70, 'cl000_3': 70, 'cl000_4': 70, 'cl000_5': 70}


    md_sim_list = ['1' for i in range(len(cl_frames))]

    base_form = '_'.join(['GAPDH_WT_NAD', ligand, 'docked'])
    #base_form = '_'.join(['GAPDH_WT_NAD_rem3PG', ligand, 'docked'])
    base_dir = '/'.join(['/home/marta', enzyme])
    scripts_dir = '/home/marta/scripts/mmpbsa_scripts'
    cpp_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']
    traj_format_list = ['crd' for i in range(len(cl_frames))]
    traj_format_i = 0
    traj_type = '10'

    for i, (cl, n_frames) in enumerate(cl_frames.items()):
        traj_format = traj_format_list[i]
        md_traj_name = ''.join([base_form, '_', cl, '_', md_sim_list[i], '_', traj_type, '.', traj_format]) if md_sim_list[i] else \
                                ''.join([base_form, '_', cl, '_', traj_type, '.', traj_format])

        md_traj_dir = '/'.join([base_dir, '3_MD_post_dock2016', form_folder, ligand_folder, cl])
        md_traj = '/'.join([md_traj_dir, md_traj_name])

        prmtop_file = '.'.join([base_form, 'prmtop'])
        leaprc_base = '/'.join([scripts_dir, 'base_files/base_leaprc'])
        mmpbsa_in_file = '/'.join([scripts_dir, 'base_files/mmpbsa_pb.in'])

        mmpbsa_base_dir = '/'.join([base_dir, '4_MMPBSA2016', form_folder, ligand_folder])

        if per_frame:
            mmpbsa_dir = '/'.join([mmpbsa_base_dir,  cl, 'per_frame'])
            if not os.path.exists(mmpbsa_dir):
                os.makedirs(mmpbsa_dir)
            prep_get_frame_from_traj(mmpbsa_dir, md_traj_dir, prmtop_file, cpp_names)

            for frame in range(1, n_frames+1):
                os.chdir(mmpbsa_dir)
                frame_dir = ''.join([mmpbsa_dir, '/frame_', str(frame)])

                # generate cpptraj files
                ligand = 'NAD,DPG'
                #ligand = 'NAD,G3P'
                prep_cpptraj_files(mmpbsa_dir, md_traj_name, ligand, cpp_names, frame)

                # generate pdb files for complex_solv/vac receptor/ligand_vac
                os.chdir(mmpbsa_dir)
                os.chmod('./run_cpptraj_mmpbsa.sh', 0755)
                subprocess.call('./run_cpptraj_mmpbsa.sh', shell=True)

                # prep the complex for mmpbsa
                #ligands = ['G3P', 'NAD']
                ligands = ['13DPG', 'NAD']

                #ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '1_G3P']),
                #                      '/'.join([base_dir, '0_Ligands_param', '3_NAD'])]
                ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '2_13DPG']),
                                     '/'.join([base_dir, '0_Ligands_param', '5_NADH'])]
                #ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '5_NADH'])]
                prep_complex_pdbs(mmpbsa_base_dir, mmpbsa_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params_list, ligands, frame)
                #prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params, ligand, frame)

                # generate an mmpbsa.in file for the current frame
                mmpbsa_dict = pb_perframe
                mmpbsa_dict['filename'] = frame_dir + '/mmpbsa_pb.in'
                mmpbsa_dict['start_frame'] = frame
                mmpbsa_dict['end_frame'] = frame

                template_folder = '/'.join([scripts_dir, 'base_files'])
                mmpbsa_per_frame_template = 'mmpbsa_perframe.jinja'
                prep_mmpbsa_in(template_folder, mmpbsa_per_frame_template, mmpbsa_dict)

                # prepare script to run mmpbsa for the current frame
                prep_mmpbsa_run(mmpbsa_dir, mmpbsa_in_file, md_traj, frame)

                # run mmpbsa for the current frame
                os.chdir(frame_dir)
                os.chmod('run_mmpbsa.sh', 0755)
                subprocess.call('./run_mmpbsa.sh', shell=True)

                # remove all files except for the results file
                try:
                    os.remove('_*')
                except OSError:
                    pass

                os.remove('mmpbsa_pb.in')
                os.remove('run_mmpbsa.sh')
                rmtree('0_complex_prep')
                # mmpbsa takes about 1min
                #time.sleep(15)
        """
        else:

            mmpbsa_dir = '/'.join([mmpbsa_base_dir,  cl])
            os.chdir(scripts_dir)

            if discard_n_frames:
                print md_traj
                print md_traj_dir
                md_traj_name_out = ''.join([md_traj_name[:-4], '_discard7f.crd'])
                discard_first_n_frames_from_traj(scripts_dir, md_traj_dir, prmtop_file, md_traj_name, md_traj_name_out, discard_n_frames)
                md_traj_name = md_traj_name_out
                md_traj = '/'.join([md_traj_dir, md_traj_name])
                print md_traj
                print '---'


            #ligand = 'PEP'
            prep_cpptraj_files(scripts_dir, md_traj_name, ligand, cpp_names)
            prep_get_frame_from_traj(scripts_dir, md_traj_dir, prmtop_file, cpp_names)  # remove?


            os.chmod('./run_cpptraj_mmpbsa.sh', 0755)
            subprocess.call('./run_cpptraj_mmpbsa.sh', shell=True)

            # if reorder_atoms:
            #
            #     pdb_base = '/'.join([base_dir, '4_MMPBSA', form_folder, ligand_folder, cl, '0_complex_prep'])
            #     off_file = '/'.join([base_dir, '0_Ligands_param', ligand_folder, ligand])
            #     n_atoms = 15
            #
            #     for cpp_name in cpp_names:
            #         pdb_file = '/'.join([pdb_base, cpp_name])
            #         reorder_ligand_atoms(pdb_file, off_file, ligand, n_atoms)


            #ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '1_G3P'])],
            #                          '/'.join([base_dir, '0_Ligands_param', '3_NAD']),
            #                          '/'.join([base_dir, '0_Ligands_param', '8_H2PO4_opt'])]

            ligands = [ligand]
            ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '3_NAD'])]
            prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params_list, ligands)

            #ligand_params = '/'.join([base_dir, '0_Ligands_param', '2_F6P'])
            #prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params, ligand)
            prep_mmpbsa_run(mmpbsa_dir, mmpbsa_in_file, md_traj)
            os.chdir(mmpbsa_dir)

            # run mmpbsa for the current frame
            os.chmod('./run_mmpbsa.sh', 0755)
            subprocess.call('nohup ./run_mmpbsa.sh &', shell=True)
        """

        traj_format_i += 1


def mmpbsa_per_frame_talb():

    per_frame = True
    reorder_atoms = False
    discard_n_frames = 7

    enzyme = 'TalB'
    form = 'HALO'
    form_folder = '_'.join(['2', form])
    ligand = 'F6P'
    ligand_folder = '_'.join(['2', ligand])

    cl_frames = {'cl201': 70, 'cl201_2': 70, 'cl201_3': 70, 'cl201_4': 70, 'cl201_5': 70,
                 'cl311': 70, 'cl311_2': 70, 'cl311_3': 70, 'cl311_4': 70, 'cl311_5': 70,
                 'cl301': 70, 'cl301_2': 70, 'cl301_3': 70, 'cl301_4': 70, 'cl301_5': 70}

    # md_sim_list = ['1']
    md_sim_list = ['1' for i in range(len(cl_frames))]

    #base_form = '_'.join(['TALB_WT_APO', ligand, 'docked'])
    base_form = '_'.join(['TALB_WT_HALO_S7P_remS7P', ligand, 'docked'])

    base_dir = '/'.join(['/home/marta', enzyme])
    scripts_dir = '/home/marta/scripts/mmpbsa_scripts'
    cpp_names = ['complex_solv', 'complex_vac', 'ligand_vac', 'receptor_vac']
    traj_format_list = ['crd' for i in range(len(cl_frames))]
    traj_format_i = 0
    traj_type = '10'

    for i, (cl, n_frames) in enumerate(cl_frames.items()):
        traj_format = traj_format_list[i]
        md_traj_name = ''.join([base_form, '_', cl, '_', md_sim_list[i], '_', traj_type, '.', traj_format]) if md_sim_list[i] else \
                                ''.join([base_form, '_', cl, '_', traj_type, '.', traj_format])

        md_traj_dir = '/'.join([base_dir, '3_MD_post_dock2016', form_folder, ligand_folder, cl])
        md_traj = '/'.join([md_traj_dir, md_traj_name])

        prmtop_file = '.'.join([base_form, 'prmtop'])
        leaprc_base = '/'.join([scripts_dir, 'base_files/base_leaprc'])
        mmpbsa_in_file = '/'.join([scripts_dir, 'base_files/mmpbsa_pb.in'])

        mmpbsa_base_dir = '/'.join([base_dir, '4_MMPBSA2016', form_folder, ligand_folder])

        if per_frame:
            mmpbsa_dir = '/'.join([mmpbsa_base_dir,  cl, 'per_frame'])
            if not os.path.exists(mmpbsa_dir):
                os.makedirs(mmpbsa_dir)
            prep_get_frame_from_traj(mmpbsa_dir, md_traj_dir, prmtop_file, cpp_names)

            for frame in range(1, n_frames+1):
                os.chdir(mmpbsa_dir)
                frame_dir = ''.join([mmpbsa_dir, '/frame_', str(frame)])

                # generate cpptraj files
                ligand = 'F6P'
                prep_cpptraj_files(mmpbsa_dir, md_traj_name, ligand, cpp_names, frame)

                # generate pdb files for complex_solv/vac receptor/ligand_vac
                os.chdir(mmpbsa_dir)
                os.chmod('./run_cpptraj_mmpbsa.sh', 0755)
                subprocess.call('./run_cpptraj_mmpbsa.sh', shell=True)

                # prep the complex for mmpbsa
                ligands = [ligand]
                ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '2_F6P'])]
                #ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '3_S7P_linear_noopt'])]

                prep_complex_pdbs(mmpbsa_base_dir, mmpbsa_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params_list, ligands, frame)
                #prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params, ligand, frame)

                # generate an mmpbsa.in file for the current frame
                mmpbsa_dict = pb_perframe
                mmpbsa_dict['filename'] = frame_dir + '/mmpbsa_pb.in'
                mmpbsa_dict['start_frame'] = frame
                mmpbsa_dict['end_frame'] = frame

                template_folder = '/'.join([scripts_dir, 'base_files'])
                mmpbsa_per_frame_template = 'mmpbsa_perframe.jinja'
                prep_mmpbsa_in(template_folder, mmpbsa_per_frame_template, mmpbsa_dict)

                # prepare script to run mmpbsa for the current frame
                prep_mmpbsa_run(mmpbsa_dir, mmpbsa_in_file, md_traj, frame)

                # run mmpbsa for the current frame
                os.chdir(frame_dir)
                os.chmod('run_mmpbsa.sh', 0755)
                subprocess.call('./run_mmpbsa.sh', shell=True)

                # remove all files except for the results file
                try:
                    os.remove('_*')
                except OSError:
                    pass

                os.remove('mmpbsa_pb.in')
                os.remove('run_mmpbsa.sh')
                rmtree('0_complex_prep')
                # mmpbsa takes about 1min
                #time.sleep(15)

        """else:

            mmpbsa_dir = '/'.join([mmpbsa_base_dir,  cl])
            os.chdir(scripts_dir)

            if discard_n_frames:
                print md_traj
                print md_traj_dir
                md_traj_name_out = ''.join([md_traj_name[:-4], '_discard7f.crd'])
                discard_first_n_frames_from_traj(scripts_dir, md_traj_dir, prmtop_file, md_traj_name, md_traj_name_out, discard_n_frames)
                md_traj_name = md_traj_name_out
                md_traj = '/'.join([md_traj_dir, md_traj_name])
                print md_traj
                print '---'


            #ligand = 'PEP'
            prep_cpptraj_files(scripts_dir, md_traj_name, ligand, cpp_names)
            prep_get_frame_from_traj(scripts_dir, md_traj_dir, prmtop_file, cpp_names)  # remove?


            os.chmod('./run_cpptraj_mmpbsa.sh', 0755)
            subprocess.call('./run_cpptraj_mmpbsa.sh', shell=True)

            # if reorder_atoms:
            #
            #     pdb_base = '/'.join([base_dir, '4_MMPBSA', form_folder, ligand_folder, cl, '0_complex_prep'])
            #     off_file = '/'.join([base_dir, '0_Ligands_param', ligand_folder, ligand])
            #     n_atoms = 15
            #
            #     for cpp_name in cpp_names:
            #         pdb_file = '/'.join([pdb_base, cpp_name])
            #         reorder_ligand_atoms(pdb_file, off_file, ligand, n_atoms)


            #ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '1_G3P'])],
            #                          '/'.join([base_dir, '0_Ligands_param', '3_NAD']),
            #                          '/'.join([base_dir, '0_Ligands_param', '8_H2PO4_opt'])]

            ligands = [ligand]
            ligand_params_list = ['/'.join([base_dir, '0_Ligands_param', '2_F6P'])]
            prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params_list, ligands)

            #ligand_params = '/'.join([base_dir, '0_Ligands_param', '2_F6P'])
            #prep_complex_pdbs(mmpbsa_base_dir, md_traj_dir, cl, cpp_names, leaprc_base, ligand_params, ligand)
            prep_mmpbsa_run(mmpbsa_dir, mmpbsa_in_file, md_traj)
            os.chdir(mmpbsa_dir)

            # run mmpbsa for the current frame
            os.chmod('./run_mmpbsa.sh', 0755)
            subprocess.call('nohup ./run_mmpbsa.sh &', shell=True)"""

        traj_format_i += 1


if __name__ == '__main__':
    mmpbsa_per_frame_eno_ab()
