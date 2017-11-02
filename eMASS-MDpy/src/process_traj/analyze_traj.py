import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from MDAnalysis import *
from MDAnalysis.analysis.distances import dist


def _add_to_dic(dic, res):
    for i in range(0, len(res), 2):
        if res[i] not in dic:
            dic[res[i]] = [res[i+1]]
        else:
            dic[res[i]].append(res[i+1])


def _filter_dic(dic, n_entries):
    for key in dic.keys():
        if len(dic[key]) < n_entries:
            del dic[key]


def _plot_out_data(data_dic, x_data, key, file_out):

    plt.scatter(x=x_data, y=data_dic[key])

    plt.xlabel('TIME / PS')
    plt.ylabel(key)
    plt.tight_layout()

    plt.savefig(''.join([file_out, '_', key, '.png']))
    plt.close()


def process_out_files(md_dir, file_in, file_out):
    """
        Processes the out files from the MD equilibration and production runs and plots all variables.

    :param md_dir: directory where the out files are
    :param file_in: absolute path to the file to be processed
    :param file_out: base name for the files to be written
    :return: None
    """

    os.chdir(md_dir)
    if not os.path.exists('data'):
        os.makedirs('data')

    data_dic = {}
    fluctuations_dic = {}
    regex_key_value = '(1-4 \S+|\S+)\s*=\s*(\-?\d+.\d+)\s*'
    regex_line = re.compile(regex_key_value*3)

    flag_averages = False
    flag_rms = False

    # process out file
    with open(file_in, 'r') as f_in:

        line = f_in.readline()
        while line:
            res = regex_line.findall(line)

            if re.findall('A V E R A G E S', line):
                flag_averages = True
            elif re.findall('R M S  F L U C T U A T I O N S', line):
                flag_rms = True

            if res and flag_averages:
                _add_to_dic(fluctuations_dic, res[0])

            elif res and flag_rms:
                _add_to_dic(fluctuations_dic, res[0])

            elif res:
                _add_to_dic(data_dic, res[0])

            line = f_in.readline()

    # process dictionaries resulting from file processing

    _filter_dic(data_dic, len(data_dic['TIME(PS)']))

    if fluctuations_dic:
        del fluctuations_dic['NSTEP']
        del fluctuations_dic['TIME(PS)']
        _filter_dic(fluctuations_dic, 2)


        # process and plot fluctuation data
        fluctuations_df = pd.DataFrame.from_dict(fluctuations_dic, orient='index')
        fluctuations_df.columns = ['average', 'rms fluctuations']
        fluctuations_df['average'] = fluctuations_df['average'].apply(float)
        fluctuations_df['rms fluctuations'] = fluctuations_df['rms fluctuations'].apply(float)

        fluctuations_df['rms relative fluctuation'] = abs(fluctuations_df['rms fluctuations'] / fluctuations_df['average'] * 100)
        fluctuations_df['rms relative fluctuation'].plot.bar()
        plt.tight_layout()
        plt.savefig(''.join([file_out, '_fluctuations.png']))
        plt.close()

    # plot variables data
    x_data = data_dic['TIME(PS)']
    for key in data_dic:
        if key not in {'TIME(PS)', 'NSTEP'}:
            _plot_out_data(data_dic, x_data, key, file_out)

    return data_dic, fluctuations_dic


def process_out_files_min(md_dir, file_in, file_out):
    """
        Processes the out file from the minimization step and plots all variables.

    :param md_dir: directory where the out files are
    :param file_in: absolute path to the file to be processed
    :param file_out: base name for the files to be written
    :return: None
    """
    os.chdir(md_dir)
    if not os.path.exists('data'):
        os.makedirs('data')

    data_dic = {}
    fluctuations_dic = {}
    regex_key_value = '(1-4 \S+|\S+)\s*=\s*(\-?\d+.\d+)\s*'
    regex_line = re.compile(regex_key_value*3)

    # process out file
    with open(file_in, 'r') as f_in:

        line = f_in.readline()
        while line:
            res = regex_line.findall(line)
            if res:
                _add_to_dic(data_dic, res[0])

            line = f_in.readline()

    _filter_dic(data_dic, 2)

    # plot variables data
    x_data = range(len(data_dic[data_dic.keys()[0]]))
    for key in data_dic:
        _plot_out_data(data_dic, x_data, key, file_out)

    return data_dic, fluctuations_dic


def _format_atom_dist(atom_dist):
    atom_dist = str('\t'.join(str(atom_dist).split()))
    atom_dist = re.sub('\[(.*)\]', r'\1', atom_dist).strip()

    return atom_dist


def _plot_atom_dist(file_out):
    dist_df = pd.read_csv('.'.join([file_out, 'dat']), sep='\t')
    del dist_df['frame']

    dist_df.plot()
    plt.xlabel('frame number')
    plt.ylabel(r'pairwise atom distance ($\AA$)')
    plt.savefig('.'.join([file_out, 'png']))
    plt.close()


def dist_backbone_ligand(traj, prmtop, ligand, file_out):
    """
        Given an MD trajectory and a ligand it finds the distance between the carbons in the ligand and the nearby
        alpha carbons from the backbone.
        The distances are saved in a TSV file and plotted.

    :param traj: path to MD trajectory
    :param prmtop: path to prmtop file
    :param ligand: ligand name as in PDB
    :param file_out: path to output file WITHOUT FILE EXTENSION
    :return atom_pairs: a list of tuples containing the atom numbers whose distances were calculated
    """

    u = Universe(prmtop, traj)
    ligand_carbons = u.select_atoms(' '.join(['resname', ligand, 'and name C*']))

    # make sure there are at least 5 atoms to calculate distances from them to the ligand carbons
    backbone_carbons = u.select_atoms(' '.join(['backbone and name CA and around 5 resname', ligand]))
    n = 5
    while len(backbone_carbons) < 5:
        backbone_carbons = u.select_atoms(' '.join(['backbone and name CA and around', str(n), 'resname', ligand]))
        n += 1

    n_atoms_comp = min(len(ligand_carbons), len(backbone_carbons))
    ligand_carbons = ligand_carbons[:n_atoms_comp]
    backbone_carbons = backbone_carbons[:n_atoms_comp]
    atom_pairs = [(ligand_carbons[i].number, backbone_carbons[i].number) for i in range(n_atoms_comp)]

    # get atom distance on frame 0
    atom_dist_ref = dist(ligand_carbons, backbone_carbons)[2]

    with open('.'.join([file_out, 'dat']), 'w') as f_out, open('_'.join([file_out, 'relative.dat']), 'w') as f_out_rel:
        # write file header
        f_out.write(''.join(['frame\t', ''.join([str(pair) + '\t' for pair in atom_pairs]).strip(), '\n']))
        f_out_rel.write(''.join(['frame\t', ''.join([str(pair) + '\t' for pair in atom_pairs]).strip(), '\n']))

        for ts in u.trajectory:
            atom_dist = dist(ligand_carbons, backbone_carbons)

            # format distance to tab separated values
            atom_dist_abs = _format_atom_dist(atom_dist[2])
            atom_dist_rel = _format_atom_dist(atom_dist[2] - atom_dist_ref)

            f_out.write(''.join([str(ts.frame), '\t', atom_dist_abs, '\n']))
            f_out_rel.write(''.join([str(ts.frame), '\t', atom_dist_rel, '\n']))

    _plot_atom_dist(file_out)
    _plot_atom_dist('_'.join([file_out, 'relative']))

    return atom_pairs


def get_atoms_dist(traj, prmtop, atom_pairs, file_out):
    """
        From a trajectory get the distances between two atoms specified in a list, atom_pairs, for each frame of the
        trajectory.

    :param traj: path to MD trajectory
    :param prmtop: path to prmtop file
    :param atom_pairs: a list of tuples, where each tuple has the atoms indices for which we want to know the distance
    :param file_out: path to output file WITHOUT FILE EXTENSION
    :return None
    """

    u = Universe(prmtop, traj)
    atom_group_1 = [u.select_atoms(''.join(['bynum ', str(pair[0])])) for pair in atom_pairs]
    atom_group_2 = [u.select_atoms(''.join(['bynum ', str(pair[1])])) for pair in atom_pairs]

    # get atom distance on frame 0
    atom_dist_ref = np.array([dist(atom1, atom2)[2][0] for atom1, atom2 in zip(atom_group_1, atom_group_2)])

    with open('.'.join([file_out, 'dat']), 'w') as f_out, open('_'.join([file_out, 'relative.dat']), 'w') as f_out_rel:
        # write file header
        f_out.write(''.join(['frame\t', ''.join([str(pair) + '\t' for pair in atom_pairs]).strip(), '\n']))
        f_out_rel.write(''.join(['frame\t', ''.join([str(pair) + '\t' for pair in atom_pairs]).strip(), '\n']))

        for ts in u.trajectory:
            atom_dist = [dist(atom1, atom2)[2][0] for atom1, atom2 in zip(atom_group_1, atom_group_2)]

            # format distance to tab separated values
            atom_dist_abs = np.array(atom_dist)
            atom_dist_rel = atom_dist_abs - atom_dist_ref

            f_out.write(''.join([str(ts.frame), '\t', _format_atom_dist(atom_dist_abs), '\n']))
            f_out_rel.write(''.join([str(ts.frame), '\t', _format_atom_dist(atom_dist_rel), '\n']))

    _plot_atom_dist(file_out)
    _plot_atom_dist('_'.join([file_out, 'relative']))


# TODO: add representations
def write_vmd_script(traj, prmtop, file_out, ligand=None, atom_pairs=None, binding_residues=None):
    """
        Writes a tcl script to load the prmtop file and the respective MD trajectory/ies.
        It assumes the MD trajectory file extension indicates its format and is 3 characters long.
        If atoms pairs are provided, the labels will also be drawn.

    :param traj: path to an MD trajectory, or a list of paths to several MD trajectories
    :param prmtop: path to the prmtop file
    :param file_out: path to the tcl script (with file extension this time)
    :param ligand: ligand name, e.g. G6P
    :param atom_pairs: a list of tuples containing the atom numbers whose distances were calculated
    :return: None
    """

    with open(file_out, 'w') as f_out:
        f_out.write(''.join(['mol new {', prmtop, '} type {parm7} first 0 last -1 step 1 waitfor 1\n']))
        if type(traj) == list and len(traj) > 1:
            for i in range(len(traj)):
                f_out.write(''.join(['mol addfile {', traj[i], '} type {', traj[i][-3:len(traj[i])], '} first 0 last -1 step 1 waitfor 1\n']))
        else:
            f_out.write(''.join(['mol addfile {', traj, '} type {', traj[-3:len(traj)], '} first 0 last -1 step 1 waitfor 1 0\n']))
        f_out.write('animate style Once\n')

        if atom_pairs:
            for atom1, atom2 in atom_pairs:
                f_out.write(''.join(['label add Bonds 0/', str(atom1), ' 0/', str(atom2), '\n']))

        # set protein representations
        f_out.write('mol addrep 0\n')
        f_out.write('mol modselect 1 0 protein\n')
        f_out.write('mol modstyle 1 0 NewCartoon 0.300000 10.000000 4.100000 0\n')
        f_out.write('mol modcolor Name\n')

        if ligand:
            # set ligand representations
            f_out.write('mol addrep 0\n')
            f_out.write(''.join(['mol modselect 2 0 resname ', ligand, '\n']))
            f_out.write('mol modstyle 2 0 CPK 1.000000 0.300000 12.000000 12.000000\n')
            f_out.write('mol modcolor Name\n')

        if binding_residues:
            rep = '3' if ligand else '2'
            f_out.write('mol addrep 0\n')
            f_out.write(''.join(['mol modselect ', rep, ' 0 residue ', ' '.join(binding_residues), '\n']))
            f_out.write(''.join(['mol modstyle ', rep, ' 0 VDW 1.0 12.000000\n']))
            f_out.write(''.join(['mol modcolor ', rep,' 0 ColorID 11\n']))

        # hide original representation
        f_out.write('mol showrep 0 0 0\n')




# TODO: get gofr and ingofr, plot it - there is an example on MDAnalysis, leave for later

# TODO: calculate radius of gyration to ensure protein didn't open up or something? leave for later

