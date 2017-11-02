import re

import numpy as np
from MDAnalysis import Universe
from MDAnalysis.analysis.distances import dist


def _get_frames_within_bounds(f_out, atom_pairs, atom_pair_bounds, atom_pair_dist_list, frame_dict):

    frames_within_bounds = []

    for i in frame_dict.keys():
        flag = 0
        for j in range(len(atom_pairs)):

            if not (atom_pair_bounds[j][0] < atom_pair_dist_list[j][i] < atom_pair_bounds[j][1]):
                flag = 1
        if not flag:
            f_out.write(''.join([frame_dict[i], '\n']))
            frames_within_bounds.append(frame_dict[i])

    f_out.write('\n\n')

    return frames_within_bounds


def _get_atom_pair_bounds(f_out, atom_pairs, atom_pair_dist_list):

    atom_pair_bounds = []

    for i, pair in enumerate(atom_pairs):
        mean = atom_pair_dist_list[i].mean()
        std = atom_pair_dist_list[i].std()
        atom_pair_bounds.append((mean - std, mean + std))

        f_out.write(''.join(['atom pair: ', str(pair), '\n']))
        f_out.write(''.join(['mean: ', str(mean), '\n']))
        f_out.write(''.join(['std: ', str(std), '\n']))
        f_out.write(''.join(['bounds: ', str(atom_pair_bounds[i]), '\n\n']))

    f_out.write('\n')

    return atom_pair_bounds


def get_rep_frame(dock_base_folder, clusters, atom_pairs, ligand, file_out):

    frames_within_bounds_dic = {}
    with open(''.join([dock_base_folder, '/', file_out, '.dat']), 'w') as f_out:
        for cl in clusters:

            with open(''.join([dock_base_folder, '/', cl, '_frames.txt']), 'r') as f_in:
                frames = f_in.read().split('\n')[:-1]  # last element is empty

            atom_pair_dist_list = [np.zeros(len(frames)) for i in range(3)]
            frame_dict = {}

            for i, frame in enumerate(frames):
                frame_dict[i] = re.findall('-(\d+)_', frame)[0]
                u = Universe(''.join([dock_base_folder, '/', frame]))

                for j, pair in enumerate(atom_pairs):
                    ligand_atom = u.select_atoms(' '.join(['resname', ligand, 'and name ', pair[0]]))
                    enzyme_atom = u.select_atoms(' '.join(['resnum', str(pair[1][0]), 'and name ', pair[1][1]]))

                    if len(ligand_atom) > 1 or len(enzyme_atom) > 1:
                        print 'More than one atom was found for the given selection - might need to be manually checked.'
                        print ligand_atom
                        print enzyme_atom
                        return

                    atom_pair_dist_list[j][i] = dist(ligand_atom, enzyme_atom)[2][0]

            f_out.write(''.join(['cluster: ', cl, '\n\n']))
            atom_pair_bounds = _get_atom_pair_bounds(f_out, atom_pairs, atom_pair_dist_list)

            f_out.write('frames within bounds:\n')
            frames_within_bounds_dic[cl] = _get_frames_within_bounds(f_out, atom_pairs, atom_pair_bounds, atom_pair_dist_list, frame_dict)

    return frames_within_bounds_dic


def _set_representations(f_cl, mol_i, ligand_resname, binding_residues=None):

    ligand = re.findall("(\w{3})$", ligand_resname)[0]

    f_cl.write('mol addrep %s\n' % mol_i)
    f_cl.write('mol modselect 1 %s protein\n' % mol_i)
    f_cl.write('mol modstyle 1 %s NewCartoon 0.300000 10.000000 4.100000 0\n' % mol_i)
    f_cl.write('mol modcolor Name\n\n')

    f_cl.write('mol addrep %s\n' % mol_i)
    f_cl.write(''.join(['mol modselect 2 %s resname ' % mol_i, ligand, '\n']))
    f_cl.write('mol modstyle 2 %s CPK 1.000000 0.300000 12.000000 12.000000\n' % mol_i)
    f_cl.write('mol modcolor 2 %s ColorID %s\n' % (mol_i, mol_i))
    f_cl.write('mol drawframes %s 2 {0:1000}\n' % mol_i)

    f_cl.write('mol showrep %s 0 0\n\n\n' % mol_i)


def get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_resname, file_out):

    for cl in frames_within_bounds_dic.keys():

        with open(''.join([dock_base_folder, '/', file_out, '_', cl, '.tcl']), 'w') as f_out:
            f_out.write(''.join(['mol new {', dock_base_folder, '/dock_prep/',
                                protein, '-%.3d' % int(frames_within_bounds_dic[cl][0]), '_receptor_', ligand,
                                '_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1\n']))

            for frame_i in range(1, len(frames_within_bounds_dic[cl])):
                f_out.write(''.join(['mol addfile {', dock_base_folder, '/dock_prep/',
                            protein, '-%.3d' % int(frames_within_bounds_dic[cl][frame_i]), '_receptor_', ligand,
                            '_docked.pdb} type {pdb} first 0 last -1 step 10 waitfor 1 0 }\n' ]))

            _set_representations(f_out, 0, ligand_resname)