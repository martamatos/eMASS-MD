from collections import OrderedDict
from collections import defaultdict

from Bio.PDB import *
from matplotlib import *

use('Agg')
from sklearn.cluster import MeanShift, estimate_bandwidth
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


# no GUI


# define some seaborn plot settings
sns.set_style('white')
sns.set_context('talk')
colors = sns.color_palette('bright')


def _get_rep_interactions(protein, protein_pdb, chain_str, ligand_resname, ligand_resnum, atom_pairs):

    parser = PDBParser()
    structure = parser.get_structure(protein, protein_pdb)

    # 0 takes the first "model" - with x-ray structure it is always 0
    model = structure[0]

    # make sure chain is the correct one!
    chain = model[chain_str]
    print atom_pairs
    print chain
    print chain[atom_pairs[0][1][0]]
    residues = [chain[atom_pairs[0][1][0]][atom_pairs[0][1][1]],
                chain[atom_pairs[1][1][0]][atom_pairs[1][1][1]],
                chain[atom_pairs[2][1][0]][atom_pairs[2][1][1]]]

    ligand_atoms = [atom_pairs[0][0], atom_pairs[1][0], atom_pairs[2][0]]

    rep_interactions = OrderedDict()
    rep_interactions[ligand_atoms[0]] = (residues[0].get_full_id(), residues[0].get_coord())
    rep_interactions[ligand_atoms[1]] = (residues[1].get_full_id(), residues[1].get_coord())
    rep_interactions[ligand_atoms[2]] = (residues[2].get_full_id(), residues[2].get_coord())

    return rep_interactions, chain, ligand_atoms


def _log_atom_interactions(f_log, rep_interactions, rep_interactions_ref, chain, map_ligand_atoms=None):

    for atom, interaction in rep_interactions.iteritems():
        f_log.write(''.join(['LIGAND ATOM: ', atom, '\n']))
        f_log.write(''.join(['CHAIN RESIDUE: ', str(interaction[0][3][1]), '\n']))
        f_log.write(''.join(['ATOM IN CHAIN RESIDUE: ', interaction[0][4][0], '\n']))
        f_log.write(''.join(['CHAIN RESIDUE:', str(chain[interaction[0][3][1]]), '\n']))
        f_log.write(''.join(['CHAIN RESIDUE + ATOM:', str(interaction), '\n\n']))

    for atom, interaction in rep_interactions_ref.iteritems():
        f_log.write(''.join(['REF LIGAND ATOM:', atom, ' map: ', map_ligand_atoms[atom] if map_ligand_atoms else atom, '\n']))
        f_log.write(''.join(['REF CHAIN RESIDUE: ', str(interaction[0][3][1]), '\n']))
        f_log.write(''.join(['REF ATOM IN CHAIN RESIDUE: ', interaction[0][4][0], '\n']))
        f_log.write(''.join(['REF CHAIN RESIDUE:', str(chain[interaction[0][3][1]]), '\n']))
        f_log.write(''.join(['REF CHAIN RESIDUE + ATOM:', str(interaction), '\n\n']))


def _get_atom_distances(docked, idx, i, protein, chain_str, ligand, ligand_resname, ligand_resnum, rep_interactions,
                        distances_dict, non_docked_frames, bond_cutoff, frames_to_discard=None, map_ligand_atoms=None):

    parser = PDBParser()
    try:
        structure = parser.get_structure('receptor_and_docked_ligand', docked)
    except:
        non_docked_frames.append(i)
        return idx, distances_dict, non_docked_frames, frames_to_discard

    model = structure[0]
    chain = model[chain_str]

    for ligand_atom, protein_interaction in rep_interactions.iteritems():
        #for residue in chain:
        #    print residue
        # first get atom object of the ligand atom
        docked_dnc_atom_object = chain[(ligand_resname, ligand_resnum, ' ')][ligand_atom]

        # then get the interaction information for that ligand atom
        interaction = protein_interaction[0]

        # then get the atom object of the residue
        interaction_atom_object = chain[interaction[3][1]][interaction[4][0]]

        # then calculate the distance
        distance = docked_dnc_atom_object - interaction_atom_object

        # if the distance between the atoms in the given atom pair is higher than a given threshold (12), filter it out
        if distance > bond_cutoff:
            frames_to_discard.append(i)

        # frame number, DNC atom, interaction atom, distance
        distances_dict[idx]['frame'] = i
        # if condition needed for when we're calculating the atom distances in the ref structure
        distances_dict[idx]['ligand_atom'] = ligand_atom if not map_ligand_atoms else map_ligand_atoms[ligand_atom]
        distances_dict[idx]['residue_atom'] = interaction[4][0]
        distances_dict[idx]['distance'] = distance

        idx += 1

    return idx, distances_dict, non_docked_frames, frames_to_discard


def _get_clusters_meanshift(dock_base_folder, data_df, atom_pair_i, param=0.4):

    new_data_df = data_df.copy()
    new_data_df['cluster'] = np.nan

    x = data_df['distance'].tolist()
    X = np.array(zip(x, np.zeros(len(x))), dtype=np.float)

    bandwidth = estimate_bandwidth(X, quantile=param)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(X)
    labels = ms.labels_
    labels_unique = np.unique(labels)

    # order clusters
    n_clusters_ = len(labels_unique)
    ordered_clusters = []
    for k in range(n_clusters_):
        my_members = labels == k
        ordered_clusters.append((min(X[my_members, 0]), k))

    ordered_clusters = sorted(ordered_clusters, key=lambda var: var[0])
    map_ordered_clusters = {item: i for i, (ind, item) in enumerate(ordered_clusters)}

    for k in range(n_clusters_):
        my_members = labels == k
        for xxx in X[my_members, 0]:
            new_data_df.loc[new_data_df[new_data_df['distance'] == xxx].index, 'cluster'] = map_ordered_clusters[k]

    new_data_df.sort_values(by='cluster', inplace=True)
    new_data_df.to_csv(''.join([dock_base_folder, '/atom_pair_', str(atom_pair_i)]), sep='\t')

    return new_data_df


def _get_cluster_limits(dock_base_folder, atom_pair, cl_i, f_log):

    f_log.write(''.join(['atom_pair', str(cl_i), '_cl.head\n', str(atom_pair.head(15)), '\n\n']))
    f_log.write(''.join(['atom_pair', str(cl_i), '_cl.cluster.unique()\n', str(atom_pair.cluster.unique()), '\n\n']))

    limits = []
    for cluster in atom_pair.cluster.unique():
        bottom = min(atom_pair[atom_pair.cluster == cluster].distance)
        top = max(atom_pair[atom_pair.cluster == cluster].distance)
        limits.append((bottom, top))

    limits = sorted(limits, key=lambda var: var[0])
    f_log.write(''.join(['limits', str(cl_i), '\n', str(limits), '\n\n']))


    atom_pair.distance.hist(bins=50, normed=1)
    i = 0
    for limit in limits:
        hist = atom_pair.distance.plot(kind='kde', grid=False,
                                       figsize=(6, 3)).axvspan(limit[0], limit[1], color=colors[i % len(colors)], alpha=0.3)
        i += 1
    fig = hist.get_figure()
    fig.savefig(''.join([dock_base_folder, '/atom', str(cl_i), '_hist.svg']), type='svg', dpi=300)
    plt.close()

    for i in range(len(limits)):
        dist_list = atom_pair[atom_pair.cluster == i].distance.values
        plt.scatter(dist_list, np.zeros_like(dist_list), color=colors[i % len(colors)])
    plt.savefig(''.join([dock_base_folder, '/atom', str(cl_i), '_1D.svg']), type='svg', dpi=300)
    plt.close()

    return limits


def _assign_frames_to_clusters_alt(distances_df, ligand_atoms, atom_pairs_cl):

    frames = distances_df.frame.unique().tolist()
    distances_df_frame_indexed = distances_df.set_index('frame')

    frames_cluster_dict = defaultdict(dict)

    for frame in frames:
        frame_cluster = ''
        for ligand_atom in distances_df_frame_indexed.loc[frame].ligand_atom.values:
            if ligand_atom == ligand_atoms[0]:
                frame_cluster += '%.0f' % atom_pairs_cl[0][atom_pairs_cl[0].frame == frame].cluster.unique()[0]

            if ligand_atom == ligand_atoms[1]:
                frame_cluster += '%.0f' % atom_pairs_cl[1][atom_pairs_cl[1].frame == frame].cluster.unique()[0]

            if ligand_atom == ligand_atoms[2]:
                frame_cluster += '%.0f' % atom_pairs_cl[2][atom_pairs_cl[2].frame == frame].cluster.unique()[0]
        frames_cluster_dict[frame]['frame_cluster'] = frame_cluster

    frames_cluster_df = pd.DataFrame.from_dict(frames_cluster_dict, orient='index')

    return frames_cluster_df


def _write_clusters(dock_base_folder, frames_cluster_df, protein, ligand, cluster_and_counts):

    for clus in cluster_and_counts.keys():
        with open(''.join([dock_base_folder, '/cl', str(clus), '_frames.txt']), 'w+') as f_cl:
            for frame in frames_cluster_df[frames_cluster_df.frame_cluster == clus].index.tolist():
                if frame != 0:
                    f_cl.write(''.join(['dock_prep/', protein, '-%.3d_receptor_' % (frame), ligand, '_docked.pdb\n']))


def _set_representations(f_cl, mol_i, ligand_resname, binding_residues=None):
    #ligand_resname = '2PG'

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

    if binding_residues:
        f_cl.write('mol addrep %s\n' % mol_i)
        f_cl.write(''.join(['mol modselect 3 %s residue ', ' '.join(binding_residues), '\n']) % mol_i)
        f_cl.write('mol modstyle 3 %s CPK 1.000000 0.300000 12.000000 12.000000\n' % mol_i)
        f_cl.write('mol modcolor 3 %s ColorID 7\n\n' % mol_i)

    f_cl.write('mol showrep %s 0 0\n\n\n' % mol_i)


def _viz_clusters_vmd(dock_base_folder, protein, ligand, ligand_resname, frames_cluster_df, cluster_and_counts,
                      ordered_clusters, ref_pdb=None, binding_residues=None):

    with open(''.join([dock_base_folder, '/viz_clusters.tcl']), 'w+') as f_cl:
        i = 0
        for clus in ordered_clusters:

            frame = frames_cluster_df[frames_cluster_df.frame_cluster == clus].index.tolist()[0]
            if ref_pdb and frame == 0:  # ref_pdb will be frame 0 and i don't want to load it as part of a cluster
                try:
                    frame = frames_cluster_df[frames_cluster_df.frame_cluster == clus].index.tolist()[1]
                except IndexError:
                    continue

            f_cl.write(''.join(['mol new {', dock_base_folder, '/dock_prep/',
                                protein, '-%.3d' % frame, '_receptor_', ligand,
                                '_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1\n']))

            f_cl.write('set lst%s' % clus)
            f_cl.write(' {')
            for x in frames_cluster_df[frames_cluster_df.frame_cluster == clus].index.tolist()[1:]:
                f_cl.write("".join(['dock_prep/', protein, '-%.3d_receptor_' % (x), ligand, '_docked.pdb ']))
            f_cl.write('}\n')
            f_cl.write('\n')
            f_cl.write('foreach i $lst%s { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 %s }\n' % (clus, i))

            _set_representations(f_cl, i, ligand_resname, binding_residues)
            i += 1

        if ref_pdb:
            f_cl.write(''.join(['mol new {', ref_pdb, '} type {pdb} first 0 last -1 step 1 waitfor 1\n']))
            _set_representations(f_cl, i, ligand_resname, binding_residues)


def cluster_docking_poses(dock_base_folder, ref_pdb, frame, protein, chain_str, n_frames_docked, ligand, ligand_resnum,
                          bond_cutoff=12, bandwidth_quantile=0.2, atom_pairs=None, ligand_resname=None, ligand_ref=None,
                          ligand_resname_ref=None, ligand_resnum_ref=None, atom_pairs_ref=None, chain_ref=None,
                          map_ligand_atoms=None,binding_residues=None, logfile=None):

    """
    Given a set of frames from a trajectory where a ligand is bound to an enzyme, cluster these frames according to
     the ligand's docking pose.

     In detail:
        1. choose 3 atom pairs (1 atom in the ligand and one in the enzyme)
        2. get the distance between the atoms in the pair
        3. for each atom pair, cluster their distances
        4. for each atom pair, determine which one is the best cluster (the one with smaller distance between the atoms)
        5. assign each frame to a cluster based on the distance between the chosen atom pairs:
            5.1 - if the given atom pair distance falls within the best cluster, it goes to cluster 1, otherwise to
                  cluster 0
            5.2 - the final cluster have the form, 101, 111, 100, etc, based on whether or not each distance between the
                  atoms in each pair belong to the best cluster or not.

    :param dock_base_folder: path to where the trajectory frames are.
    :param ref_pdb: the pdb file to take as a reference for the docking
    :param frame: the first frame of the trajectory
    :param protein: the protein name, e.g. TALB_WT_HALO_S7P_remS7P
    :param chain_str: the protein chain to be considered
    :param n_frames_docked: how many docked frames exist
    :param ligand: the ligand's name as in the system's folders
    :param ligand_resnum: the residue number of the ligand in the docked frames
    :param bond_cutoff: if the distance between the atoms in an atom_pair is greater than this, the frame will be discarded
    :param bandwidth_quantile:
    :param atom_pairs: the atom pairs to be considered, have the form
                            ('ligand_atom', (enzyme_residue_number, 'residue_atom_name')), e.g. ('O2', (131, 'NZ'))
    :param ligand_resname: the residue name of the ligand, e.g. H_G6P
    :param ligand_ref: the ligand name in the reference pdb (ref_pdb)
    :param ligand_resname_ref: the ligand residue name in the reference pdb (ref_pdb)
    :param ligand_resnum_ref: the ligand residue number in the reference pdb (ref_pdB)
    :param atom_pairs_ref: the atom pairs in the reference_pdb for which the distances will be calculated
    :param chain_ref: the enzyme chain to be considered in the reference pdb (ref_pdb)
    :param map_ligand_atoms: a dictionary mapping the ligand atoms in the docked frames to the ligand atoms in
                             reference pdb, e.g. {'O2': 'O5'}
    :param binding_residues: a list with the residue numbers of the binding residues in the enzyme
    :param logfile: path to a file where the output from clustering will be logged
    :return: a list with the names of the 3 largest clusters
    """

    if not ligand_ref:
        ligand_ref = ligand

    if not ligand_resnum_ref:
        ligand_resnum_ref = ligand_resnum

    if not ligand_resname:
        ligand_resname = '_'.join(['H', ligand])

    if not ligand_resname_ref:
        ligand_resname_ref = '_'.join(['H', ligand_ref])

    if not atom_pairs_ref:
        atom_pairs_ref = atom_pairs

    if not chain_ref:
        chain_ref = chain_str

    if not logfile:
        logfile = '/'.join([dock_base_folder, 'log_clustering.txt'])

    with open(logfile, 'w') as f_log:

        rep_interactions_ref, chain, ligand_atoms = _get_rep_interactions(protein, ref_pdb, chain_ref,
                                                                          ligand_resname_ref, ligand_resnum_ref,
                                                                          atom_pairs_ref)
        rep_interactions, chain, ligand_atoms = _get_rep_interactions(protein, frame, chain_str, ligand_resname,
                                                                      ligand_resnum, atom_pairs)

        _log_atom_interactions(f_log, rep_interactions, rep_interactions_ref, chain, map_ligand_atoms)

        # get ref atom pair distances
        idx = 0
        distances_dict = defaultdict(dict)
        non_docked_frames = []
        frames_to_discard = []
        idx, distances_dict, non_docked_frames, frames_to_discard =\
            _get_atom_distances(ref_pdb, idx, 0, protein, chain_ref, ligand_ref, ligand_resname_ref, ligand_resnum_ref,
                                rep_interactions_ref, distances_dict, non_docked_frames, bond_cutoff, frames_to_discard,
                                map_ligand_atoms=map_ligand_atoms)

        if frames_to_discard:
            distances_dict = defaultdict(dict)
            idx = 0

        # get all frames atom pair distances
        frames_to_discard = []
        for i in range(1, n_frames_docked+1):

            docked = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % i + ligand + '_docked.pdb'
            idx, distances_dict, non_docked_frames, frames_to_discard = \
                _get_atom_distances(docked, idx, i, protein, chain_str, ligand, ligand_resname, ligand_resnum,
                                    rep_interactions, distances_dict, non_docked_frames, bond_cutoff, frames_to_discard)

        frames_to_discard = list(set(frames_to_discard))  # remove duplicates

        f_log.write(''.join(['number of frames not docked:\n', str(len(non_docked_frames)), '\n\n']))
        f_log.write(''.join(['frames not docked:\n', str(non_docked_frames), '\n\n']))
        f_log.write(''.join(['number of frames to discard:\n', str(len(frames_to_discard)), '\n']))
        f_log.write(''.join(['frames to discard:\n', str(frames_to_discard), '\n\n']))


        distances_df = pd.DataFrame.from_dict(distances_dict, orient='index')
        distances_df = distances_df[['frame', 'ligand_atom', 'residue_atom', 'distance']]
        distances_df.to_csv('/'.join([dock_base_folder, 'distances']), sep='\t')

        for frame in frames_to_discard:
            distances_df = distances_df[distances_df.frame != frame]


        #  Clustering the distances between each atom pair
        # Ligand atom #1 to residue atom #1
        atom_pair1_cl = _get_clusters_meanshift(dock_base_folder, distances_df[distances_df.ligand_atom == ligand_atoms[0]],
                                                1, bandwidth_quantile)
        _get_cluster_limits(dock_base_folder, atom_pair1_cl, 1, f_log)

        # Ligand atom #2 to residue atom #2
        atom_pair2_cl = _get_clusters_meanshift(dock_base_folder, distances_df[distances_df.ligand_atom == ligand_atoms[1]],
                                                2, bandwidth_quantile)
        _get_cluster_limits(dock_base_folder, atom_pair2_cl, 2, f_log)

        # Ligand atom #3 to residue atom #3
        atom_pair3_cl = _get_clusters_meanshift(dock_base_folder, distances_df[distances_df.ligand_atom == ligand_atoms[2]],
                                                3, bandwidth_quantile)
        _get_cluster_limits(dock_base_folder, atom_pair3_cl, 3, f_log)


        # assigning frames to clusters
        atom_pairs_cl = [atom_pair1_cl, atom_pair2_cl, atom_pair3_cl]
        frames_cluster_df = _assign_frames_to_clusters_alt(distances_df, ligand_atoms, atom_pairs_cl)

        cluster_and_counts = frames_cluster_df.frame_cluster.value_counts().to_dict()
        ordered_clusters = frames_cluster_df.frame_cluster.value_counts().index.values

        f_log.write(''.join(['frames_cluster_df.head()\n', str(frames_cluster_df.head()), '\n\n']))
        f_log.write(''.join(['frames_cluster_df.frame_cluster.value_counts()/len(frames_cluster_df)\n',
                             str(frames_cluster_df.frame_cluster.value_counts()/len(frames_cluster_df)), '\n\n']))
        f_log.write(''.join(['CLUSTERS FOR MMPBSA:\n', str(ordered_clusters)]))

    # create file with frames in each cluster
    _write_clusters(dock_base_folder, frames_cluster_df, protein, ligand, cluster_and_counts)

    # visualizing clusters in VMD, write tcl files to load clusters
    _viz_clusters_vmd(dock_base_folder, protein, ligand, ligand_resname, frames_cluster_df, cluster_and_counts,
                      ordered_clusters, ref_pdb, binding_residues)

    return ['cl' + item for item in ordered_clusters[:5]]


