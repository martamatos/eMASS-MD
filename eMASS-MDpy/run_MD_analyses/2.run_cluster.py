from src.cluster_dock.cluster_dock import cluster_docking_poses
from src.cluster_dock.analyze_clusters import get_rep_frame, get_vmd_viz_rep_frames


def cluster_ENO_APO_2PG():

    dock_base_folder = '/home/mrama/Desktop/MD/ENO_1E9I/2_Dock/1_APO/1_2PG'

    chain = 'A'
    chain_ref = ' '
    protein = '1E9I_WT_APO'
    ligand = '2PG'
    ligand_ref = '2PG'
    #ligand_resname = 'H_2PG'
    ligand_resname_ref = 'H_2PG'

    ligand_resnum = 434
    ligand_resnum_ref = 433
    n_pdbs_docked = 750  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/ENO_1E9I/0_Crystal_structures/3H8A_A_WT_2PG.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C3', (208, 'CD')), ('P1', (370, 'CZ')), ('C1', (341, 'NZ'))]
    atom_pairs_ref = [('C3', (208, 'CD')), ('P', (370, 'CZ')), ('C2', (341, 'NZ'))]  # 3.5, 4.85,  2.7
    map_ligand_atoms = {'C3': 'C3', 'P': 'P1', 'C2': 'C1'}
    bandwidth_quantile = 0.15
    bond_cutoff = 15

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=None, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_ENO_APO_PEP():

    dock_base_folder = '/home/mrama/Desktop/MD/ENO_1E9I/2_Dock/1_APO/2_PEP'

    chain = 'A'
    chain_ref = ' '
    protein = '1E9I_WT_APO'
    ligand = 'PEP'
    ligand_ref = 'PEP'
    #ligand_resname = 'H_2PG'
    ligand_resname_ref = 'H_PEP'

    ligand_resnum = 434
    ligand_resnum_ref = 433
    n_pdbs_docked = 750  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/ENO_1E9I/0_Crystal_structures/3H8A_A_WT_PEP.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C3', (208, 'CD')), ('P', (370, 'CZ')), ('C2', (341, 'NZ'))]
    atom_pairs_ref = [('C3', (208, 'CD')), ('P', (370, 'CZ')), ('C2', (341, 'NZ'))]  # 3.5, 4.85, 2.7
    map_ligand_atoms = {'C3': 'C3', 'C2': 'C2', 'P': 'P'}
    bandwidth_quantile = 0.17
    bond_cutoff = 15

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=None, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)



def cluster_ENO_APO_2PG_liz():

    dock_base_folder = '/home/mrama/Desktop/MD/ENO_AB/2_Dock/1_MG/1_2PG'

    chain = ' '
    chain_ref = ' '
    protein = 'ENO_AB_WT_MG'
    ligand = '2PG'
    ligand_ref = '2PG'
    ligand_resname = ' '
    ligand_resname_ref = 'H_2PG'

    ligand_resnum = 433
    ligand_resnum_ref = 433
    n_pdbs_docked = 252  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/ENO_1E9I/0_Crystal_structures/3H8A_A_WT_2PG.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C3', (208, 'CD')), ('P', (370, 'CZ')), ('C1', (341, 'NZ'))]
    atom_pairs_ref = [('C3', (208, 'CD')), ('P', (370, 'CZ')), ('C2', (341, 'NZ'))]  # 3.5, 4.85,  2.7
    map_ligand_atoms = {'C3': 'C3', 'P': 'P', 'C2': 'C1'}
    bandwidth_quantile = 0.28
    bond_cutoff = 15

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_ENO_APO_PEP_liz():

    dock_base_folder = '/home/mrama/Desktop/MD/ENO_AB/2_Dock/1_MG/2_PEP'

    chain = ' '
    chain_ref = ' '
    protein = 'ENO_AB_WT_MG'
    ligand = 'PEP'
    ligand_ref = 'PEP'
    ligand_resname = ' '
    ligand_resname_ref = 'H_PEP'

    ligand_resnum = 433
    ligand_resnum_ref = 433
    n_pdbs_docked = 252  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/ENO_1E9I/0_Crystal_structures/3H8A_A_WT_PEP.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C3', (208, 'CD')), ('P', (370, 'CZ')), ('C1', (341, 'NZ'))]
    atom_pairs_ref = [('C3', (208, 'CD')), ('P', (370, 'CZ')), ('C2', (341, 'NZ'))]  # 3.5, 4.85,  2.7
    map_ligand_atoms = {'C3': 'C3', 'P': 'P', 'C2': 'C1'}
    bandwidth_quantile = 0.27
    bond_cutoff = 15

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)



def cluster_GAPDH_APO_G3P():

    dock_base_folder = '/home/mrama/Desktop/MD/GAPDH/2_Dock_reproduce/1_APO/1_G3P'

    chain = 'A'
    chain_ref = 'A'
    protein = 'GAPDH_WT_APO'
    ligand = 'G3P'
    ligand_ref = 'G3H'
    ligand_resname = 'H_G3H'

    ligand_resnum = 331
    ligand_resnum_ref = 350
    n_pdbs_docked = 750  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/GAPDH/0_Crystal_structures/1DC4_chainA.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('O1', (209, 'CA')), ('C1', (176, 'CE1')), ('C2', (149, 'SG'))]
    atom_pairs_ref = [('O1P', (209, 'CA')), ('C3', (176, 'CE1')), ('C2', (149, 'SG'))]
    map_ligand_atoms = {'O1P': 'O1', 'C3': 'C1', 'C2': 'C2'}
    bandwidth_quantile = 0.2
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref,ligand_resname_ref=None,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_GAPDH_APO_DPG():

    dock_base_folder = '/home/mrama/Desktop/MD/GAPDH/2_Dock_reproduce/1_APO/2_13DPG'

    chain = 'A'
    chain_ref = 'A'
    protein = 'GAPDH_WT_APO'
    ligand = '13DPG'
    ligand_ref = 'DPG'
    ligand_resname = 'H_DPG'

    ligand_resnum = 331
    ligand_resnum_ref = 331
    n_pdbs_docked = 750  # consider just counting the number of files in the directory
    ref_pdb = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('O4', (209, 'CA')), ('C1', (176, 'CE1')), ('C2', (149, 'SG'))]
    map_ligand_atoms = {'O4': 'O4', 'C2': 'C2', 'C1': 'C1'}
    bandwidth_quantile = 0.2
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=None,
                                         ligand_resnum_ref=ligand_resnum_ref,atom_pairs_ref=None, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_GAPDH_NAD_G3P():

    dock_base_folder = '/home/mrama/Desktop/MD/GAPDH/2_Dock_reproduce/2_NAD/1_G3P'

    chain = 'A'
    chain_ref = 'A'
    protein = 'GAPDH_WT_NAD'
    ligand = 'G3P'
    ligand_ref = 'G3H'
    ligand_resname = 'H_G3H'
    #ligand_resname_ref = 'G3H'

    ligand_resnum = 332
    ligand_resnum_ref = 350
    n_pdbs_docked = 750  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/GAPDH/0_Crystal_structures/1DC4_chainA.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    # GAPDH, APO, G3P
    atom_pairs = [('O1', (209, 'CA')), ('C1', (176, 'CE1')), ('C2', (149, 'SG'))]
    atom_pairs_ref = [('O1P', (209, 'CA')), ('C3', (176, 'CE1')), ('C2', (149, 'SG'))]
    map_ligand_atoms = {'O1P': 'O1', 'C3': 'C1', 'C2': 'C2'}
    bandwidth_quantile = 0.2
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=None,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_GAPDH_NAD_DPG():

    dock_base_folder = '/home/mrama/Desktop/MD/GAPDH/2_Dock_reproduce/2_NAD/2_13DPG'

    chain = 'A'
    chain_ref = 'A'
    protein = 'GAPDH_WT_NAD'
    ligand = '13DPG'
    ligand_ref = 'DPG'
    ligand_resname = 'H_DPG'

    ligand_resnum = 332
    ligand_resnum_ref = 332
    n_pdbs_docked = 750  # consider just counting the number of files in the directory
    ref_pdb = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    # GAPDH, APO, G3P
    atom_pairs = [('O4', (209, 'CA')), ('C1', (176, 'CE1')), ('C2', (149, 'SG'))]
    map_ligand_atoms = {'O4': 'O4', 'C1': 'C1', 'C2': 'C2'}
    bandwidth_quantile = 0.2
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=None,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=None, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_GAPDH_NAD_DPG_based_on_crystal():

    dock_base_folder = '/home/mrama/Desktop/MD/GAPDH/2_Dock_reproduce/2_NAD/2_13DPG'

    chain = 'A'
    chain_ref = 'A'
    protein = 'GAPDH_WT_NAD'
    ligand = '13DPG'
    ligand_ref = 'G3H'
    ligand_resname = 'H_DPG'
    ligand_resname_ref = 'H_G3H'

    ligand_resnum = 332
    ligand_resnum_ref = 350
    n_pdbs_docked = 750  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/GAPDH/0_Crystal_structures/1DC4_chainA.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    # GAPDH, APO, G3P
    atom_pairs = [('C1', (149, 'SG')), ('P2', (209, 'CA')), ('C3', (231, 'CZ'))]
    atom_pairs_ref = [('C1', (149, 'SG')), ('P', (209, 'CA')), ('C3', (231, 'CZ'))]
    # CYS149:SG/G3H350:C1 = 1.824102, GLY209:CA/G3H350:P = 4.633661, ARG231:CZ/G3H350:C3 = 4.826069
    map_ligand_atoms = {'C1': 'C1', 'P': 'P2', 'C3': 'C3'}
    bandwidth_quantile = 0.25
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    top_clusters = ['cl101', 'cl100', 'cl220']
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, 'DPG', file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, 'DPG', ligand_ref, file_out)


def cluster_GAPDH_HALO_remG3P_G3P():

    dock_base_folder = '/home/mrama/Desktop/MD/GAPDH/2_Dock_reproduce/4_NAD_G3P_remG3P/1_G3P'

    chain = 'A'
    chain_ref = 'A'
    protein = 'GAPDH_WT_NAD_rem3PG'
    ligand = 'G3P'
    ligand_ref = 'G3H'
    ligand_resname = 'H_G3H'
    #ligand_resname_ref = 'G3H'

    ligand_resnum = 332
    ligand_resnum_ref = 350
    n_pdbs_docked = 745  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/GAPDH/0_Crystal_structures/1DC4_chainA.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    # GAPDH, APO, G3P
    atom_pairs = [('O1', (209, 'CA')), ('C1', (176, 'CE1')), ('C2', (149, 'SG'))]
    atom_pairs_ref = [('O1P', (209, 'CA')), ('C3', (176, 'CE1')), ('C2', (149, 'SG'))]
    map_ligand_atoms = {'O1P': 'O1', 'C3': 'C1', 'C2': 'C2'}
    bandwidth_quantile = 0.2
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=None,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref,chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_GAPDH_HALO_remG3P_DPG():

    dock_base_folder = '/home/mrama/Desktop/MD/GAPDH/2_Dock_reproduce/4_NAD_G3P_remG3P/2_13DPG_test'

    chain = 'A'
    chain_ref = 'A'
    protein = 'GAPDH_WT_NAD_rem3PG'
    ligand = '13DPG'
    ligand_ref = 'DPG'
    ligand_resname = 'H_DPG'

    ligand_resnum = 332
    ligand_resnum_ref = 332
    n_pdbs_docked = 745  # consider just counting the number of files in the directory
    ref_pdb = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('O4', (209, 'CA')), ('C1', (176, 'CE1')), ('C2', (149, 'SG'))]
    map_ligand_atoms = {'O4': 'O4', 'C1': 'C1', 'C2': 'C2'}
    bandwidth_quantile = 0.2
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=None,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=None, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)



def cluster_GAPDH_HALO_remG3P_DPG_based_on_crystal():

    dock_base_folder = '/home/mrama/Desktop/MD/GAPDH/2_Dock_reproduce/4_NAD_G3P_remG3P/2_13DPG'

    chain = 'A'
    chain_ref = 'A'
    protein = 'GAPDH_WT_NAD_rem3PG'
    ligand = '13DPG'
    ligand_ref = 'G3H'
    ligand_resname = 'H_DPG'
    ligand_resname_ref = 'H_G3H'

    ligand_resnum = 332
    ligand_resnum_ref = 350
    n_pdbs_docked = 745  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/GAPDH/0_Crystal_structures/1DC4_chainA.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    # GAPDH, APO, G3P
    atom_pairs = [('C1', (149, 'SG')), ('P2', (209, 'CA')), ('C3', (231, 'CZ'))]
    atom_pairs_ref = [('C1', (149, 'SG')), ('P', (209, 'CA')), ('C3', (231, 'CZ'))]
    # CYS149:SG/G3H350:C1 = 1.824102, GLY209:CA/G3H350:P = 4.633661, ARG231:CZ/G3H350:C3 = 4.826069
    map_ligand_atoms = {'C1': 'C1', 'P': 'P2', 'C3': 'C3'}
    bandwidth_quantile = 0.23
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    #top_clusters = ['cl101', 'cl100', 'cl220']
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, 'DPG', file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, 'DPG', ligand_ref, file_out)


def cluster_TALB_APO_S7P():

    dock_base_folder = '/home/mrama/Desktop/MD/TALB/2_Dock_reproduce/1_APO/1_S7P'

    chain = 'A'
    chain_ref = ' '
    protein = 'TALB_WT_APO'
    ligand = 'S7P'
    ligand_ref = 'S7P'
    #ligand_resname = 'H_G6P'
    ligand_resname_ref = ' '

    ligand_resnum = 317
    ligand_resnum_ref = 317
    n_pdbs_docked = 803  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/TALB/0_Crystalline_structures/1UCW_wt_S7P_watbox.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C1', (32, 'CG2')), ('C', (131, 'CD')), ('C3', (242, 'CG2'))]
    atom_pairs_ref = [('C2', (32, 'CG2')), ('C1', (131, 'CD')), ('C4', (242, 'CG2'))]
    map_ligand_atoms = {'C2': 'C1', 'C1': 'C', 'C4': 'C3'}
    bandwidth_quantile = 0.15
    bond_cutoff = 20

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=None, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_TALB_APO_S7P_linear():

    dock_base_folder = '/home/mrama/Desktop/MD/TALB/2_Dock_reproduce/1_APO/1_S7P_linear'

    chain = 'A'
    chain_ref = ' '
    protein = 'TALB_WT_APO'
    ligand = 'S7P'
    ligand_ref = 'S7P'
    #ligand_resname = 'H_G6P'
    ligand_resname_ref = ' '

    ligand_resnum = 317
    ligand_resnum_ref = 317
    n_pdbs_docked = 803  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/TALB/0_Crystalline_structures/1UCW_wt_S7P_watbox.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C2', (32, 'CG2')), ('C3', (131, 'CD')), ('C4', (242, 'CG2'))]
    atom_pairs_ref = [('C2', (32, 'CG2')), ('C1', (131, 'CD')), ('C4', (242, 'CG2'))]
    map_ligand_atoms = {'C2': 'C2', 'C1': 'C3', 'C4': 'C4'}
    bandwidth_quantile = 0.25
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=None, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_TALB_APO_F6P():

    dock_base_folder = '/home/mrama/Desktop/MD/TalB/2_Dock_reproduce/1_APO/2_F6P'

    chain = 'A'
    chain_ref = ' '
    protein = 'TALB_WT_APO'
    ligand = 'F6P'
    ligand_ref = 'S7P'
    ligand_resname = 'H_F6P'
    ligand_resname_ref = ' '

    ligand_resnum = 317
    ligand_resnum_ref = 317
    n_pdbs_docked = 803  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/TalB/0_Crystalline_structures/1UCW_wt_S7P_watbox.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C5', (32, 'CG2')), ('C4', (131, 'CD')), ('C3', (242, 'CG2'))]
    atom_pairs_ref = [('C2', (32, 'CG2')), ('C1', (131, 'CD')), ('C4', (242, 'CG2'))]
    map_ligand_atoms = {'C2': 'C5', 'C1': 'C4', 'C4': 'C3'}
    #bandwidth_quantile = 0.23
    bandwidth_quantile = 0.22
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref,atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, 'F6P', file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_TALB_APO_F6P_linear():

    dock_base_folder = '/home/mrama/Desktop/MD/TALB/2_Dock_reproduce/1_APO/2_F6P'

    chain = 'A'
    chain_ref = 'A'
    protein = 'TALB_WT_APO'
    ligand = 'F6P'
    ligand_ref = 'F6P'
    #ligand_resname = 'H_G6P'
    #ligand_resname_ref = 'F6P'

    ligand_resnum = 317
    ligand_resnum_ref = 317
    n_pdbs_docked = 803  # consider just counting the number of files in the directory
    ref_pdb = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C5', (32, 'CG2')), ('C4', (131, 'CD')), ('C3', (242, 'CG2'))]
    map_ligand_atoms = {'C5': 'C5', 'C4': 'C4', 'C3': 'C3'}
    bandwidth_quantile = 0.2
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=None, ligand_ref=ligand_ref, ligand_resname_ref=None,
                                         ligand_resnum_ref=ligand_resnum_ref,atom_pairs_ref=atom_pairs, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_TALB_APO_F6P_linear_4S2C():

    dock_base_folder = '/home/mrama/Desktop/MD/TalB/2_Dock_reproduce/1_APO/2_F6P_4S2C'

    chain = 'A'
    chain_ref = 'A'
    protein = 'TALB_WT_APO'
    ligand = 'F6P'
    ligand_ref = 'F6R'
    ligand_resname = 'H_F6P'
    ligand_resname_ref = 'H_F6R'

    ligand_resnum = 317
    ligand_resnum_ref = 318
    n_pdbs_docked = 803  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/TalB/0_Crystalline_structures/4s2c_A.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C2', (131, 'NZ')), ('C6', (34, 'ND2')), ('P', (227, 'CZ'))]
    atom_pairs_ref = [('C2', (132, 'NZ')), ('C6', (35, 'ND2')), ('P', (228, 'CZ'))]
    map_ligand_atoms = {'C2': 'C2', 'C6': 'C6', 'P': 'P'}
    bandwidth_quantile = 0.3
    bond_cutoff = 15

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref,atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand, file_out)


def cluster_TALB_HALO_S7P():

    dock_base_folder = '/home/mrama/Desktop/MD/TALB/2_Dock_reproduce/2_HALO/1_S7P'

    chain = 'A'
    chain_ref = ' '
    protein = 'TALB_WT_HALO_S7P_remS7P'
    ligand = 'S7P'
    ligand_ref = 'S7P'
    #ligand_resname = 'H_G6P'
    ligand_resname_ref = ' '

    ligand_resnum = 317
    ligand_resnum_ref = 317
    n_pdbs_docked = 493  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/TALB/0_Crystalline_structures/1UCW_wt_S7P_watbox.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C1', (32, 'CG2')), ('C', (131, 'CD')), ('C3', (242, 'CG2'))]
    atom_pairs_ref = [('C2', (32, 'CG2')), ('C1', (131, 'CD')), ('C4', (242, 'CG2'))]
    map_ligand_atoms = {'C2': 'C1', 'C1': 'C', 'C4': 'C3'}
    #bandwidth_quantile = 0.23
    bandwidth_quantile = 0.23
    bond_cutoff = 12


    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=None, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_TALB_HALO_S7P_linear():

    dock_base_folder = '/home/mrama/Desktop/MD/TALB/2_Dock_reproduce/2_HALO/1_S7P_linear'

    chain = 'A'
    chain_ref = ' '
    protein = 'TALB_WT_HALO_S7P_remS7P'
    ligand = 'S7P'
    ligand_ref = 'S7P'
    #ligand_resname = 'H_G6P'
    ligand_resname_ref = ' '

    ligand_resnum = 317
    ligand_resnum_ref = 317
    n_pdbs_docked = 493  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/TALB/0_Crystalline_structures/1UCW_wt_S7P_watbox.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C2', (32, 'CG2')), ('C3', (131, 'CD')), ('C4', (242, 'CG2'))]
    atom_pairs_ref = [('C2', (32, 'CG2')), ('C1', (131, 'CD')), ('C4', (242, 'CG2'))]
    map_ligand_atoms = {'C2': 'C2', 'C1': 'C3', 'C4': 'C4'}
    #bandwidth_quantile = 0.23
    bandwidth_quantile = 0.2
    bond_cutoff = 12


    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=None, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref, atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_TALB_HALO_F6P():

    dock_base_folder = '/home/mrama/Desktop/MD/TalB/2_Dock_reproduce/2_HALO/2_F6P'

    chain = 'A'
    chain_ref = ' '
    protein = 'TALB_WT_HALO_S7P_remS7P'
    ligand = 'F6P'
    ligand_ref = 'S7P'
    ligand_resname = 'H_F6P'
    ligand_resname_ref = ' '

    ligand_resnum = 317
    ligand_resnum_ref = 317
    n_pdbs_docked = 492  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/TalB/0_Crystalline_structures/1UCW_wt_S7P_watbox.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C5', (32, 'CG2')), ('C4', (131, 'CD')), ('C3', (242, 'CG2'))]
    atom_pairs_ref = [('C2', (32, 'CG2')), ('C1', (131, 'CD')), ('C4', (242, 'CG2'))]
    map_ligand_atoms = {'C2': 'C5', 'C1': 'C4', 'C4': 'C3'}
    #bandwidth_quantile = 0.23
    bandwidth_quantile = 0.25
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref,atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, 'F6P', file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_TALB_HALO_F6P_linear():

    dock_base_folder = '/home/mrama/Desktop/MD/TALB/2_Dock_reproduce/2_HALO/2_F6P'

    chain = 'A'
    chain_ref = 'A'
    protein = 'TALB_WT_HALO_S7P_remS7P'
    ligand = 'F6P'
    ligand_ref = 'F6P'
    #ligand_resname = 'H_G6P'
    #ligand_resname_ref = 'F6P'

    ligand_resnum = 317
    ligand_resnum_ref = 317
    n_pdbs_docked = 493  # consider just counting the number of files in the directory
    ref_pdb = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C5', (32, 'CG2')), ('C4', (131, 'CD')), ('C3', (242, 'CG2'))]
    map_ligand_atoms = {'C5': 'C5', 'C4': 'C4', 'C3': 'C3'}
    bandwidth_quantile = 0.27
    #bandwidth_quantile = 0.2
    bond_cutoff = 12

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=None, ligand_ref=ligand_ref, ligand_resname_ref=None,
                                         ligand_resnum_ref=ligand_resnum_ref,atom_pairs_ref=atom_pairs, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand_ref, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand_ref, file_out)


def cluster_TALB_HALO_F6P_linear_4S2C():

    dock_base_folder = '/home/mrama/Desktop/MD/TalB/2_Dock_reproduce/2_HALO/2_F6P_4S2C'

    chain = 'A'
    chain_ref = 'A'
    protein = 'TALB_WT_HALO_S7P_remS7P'
    ligand = 'F6P'
    ligand_ref = 'F6R'
    ligand_resname = 'H_F6P'
    ligand_resname_ref = 'H_F6R'

    ligand_resnum = 317
    ligand_resnum_ref = 318
    n_pdbs_docked = 803  # consider just counting the number of files in the directory
    ref_pdb = '/home/mrama/Desktop/MD/TalB/0_Crystalline_structures/4s2c_A.pdb'

    frame0 = dock_base_folder + '/dock_prep/' + protein + '-%.3d_receptor_' % 1 + ligand + '_docked.pdb'

    atom_pairs = [('C2', (131, 'NZ')), ('C6', (34, 'ND2')), ('P', (227, 'CZ'))]
    atom_pairs_ref = [('C2', (132, 'NZ')), ('C6', (35, 'ND2')), ('P', (228, 'CZ'))]
    map_ligand_atoms = {'C2': 'C2', 'C6': 'C6', 'P': 'P'}
    bandwidth_quantile = 0.4
    bond_cutoff = 20

    top_clusters = cluster_docking_poses(dock_base_folder, ref_pdb, frame0, protein, chain, n_pdbs_docked, ligand, ligand_resnum,
                                         bond_cutoff=bond_cutoff, bandwidth_quantile=bandwidth_quantile, atom_pairs=atom_pairs,
                                         ligand_resname=ligand_resname, ligand_ref=ligand_ref, ligand_resname_ref=ligand_resname_ref,
                                         ligand_resnum_ref=ligand_resnum_ref,atom_pairs_ref=atom_pairs_ref, chain_ref=chain_ref,
                                         map_ligand_atoms=map_ligand_atoms)

    file_out = 'rep_frame'
    frames_within_bounds_dic = get_rep_frame(dock_base_folder, top_clusters, atom_pairs, ligand, file_out)
    get_vmd_viz_rep_frames(frames_within_bounds_dic, dock_base_folder, protein, ligand, ligand, file_out)



if __name__ == '__main__':

    #cluster_ENO_APO_2PG()
    #cluster_ENO_APO_PEP()

    # the reasoning to choose the atom pairs for DPG was to keep the pairs equivalent to the ones in
    #   the crystal structure when G3P is bound to the enzyme (chosen atoms in the ligand: core C/O)

    #cluster_GAPDH_APO_G3P()
    #cluster_GAPDH_APO_DPG()
    #cluster_GAPDH_NAD_G3P()
    #cluster_GAPDH_NAD_DPG()
    #cluster_GAPDH_HALO_remG3P_G3P()
    #cluster_GAPDH_HALO_remG3P_DPG()

    #cluster_GAPDH_NAD_DPG_based_on_crystal()

    cluster_ENO_APO_2PG_liz()

    # the reasoning to choose the atom pairs for F6P was to keep the pairs equivalent to the ones in
    #   the crystal structure when S7P is bound to the enzyme (chosen atoms in the ligand: core carbons)

    # cluster_TALB_APO_S7P()
    # cluster_TALB_APO_S7P_linear()
    # cluster_TALB_APO_F6P()
    # cluster_TALB_HALO_S7P()
    #cluster_TALB_HALO_S7P_linear()
    # cluster_TALB_HALO_F6P()

    #cluster_TALB_APO_F6P_4S2C()
    #cluster_TALB_HALO_F6P_4S2C()

    #cluster_TALB_APO_F6P()
    #cluster_TALB_HALO_F6P()
