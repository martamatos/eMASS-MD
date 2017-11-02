from MDAnalysis import Universe
from MDAnalysis.analysis.align import alignto

from src.process_traj.rmsd_MDAnalysis import calculate_traj_rmsd_fit


def align_enzymes(ref_pdb, frame_pdb, output_pdb):
    """
    Given a reference pdb and frame_pdb aligns the protein backbones.

    :param ref_pdb: path+name for reference pdb file
    :param frame_pdb: path+name for frame pdb file
    :param output_pdb: path+name for output pdb file with aligned backbones
    :return: None
    """

    ref = Universe(ref_pdb)
    frame = Universe(frame_pdb)
    alignto(frame, ref, select='backbone', mass_weighted=False)

    frame.atoms.write(output_pdb)


def generate_vmd_viz(file_out, ref_pdb, aligned_pdb, ligand):
    """
    Generates pdb file with reference pdb and pdb with aligned enzyme structures.

    :param file_out: path+name for output pdb file
    :param ref_pdb: path+name for reference pdb file
    :param aligned_pdb: path+name for pdb file with aligned structures
    :param ligand: ligand name
    :return: None
    """

    with open(file_out, 'w') as f_out:
        f_out.write(''.join(['mol new {', ref_pdb, '} type {pdb} first 0 last -1 step 1 waitfor 1\n']))
        f_out.write(''.join(['mol new {', aligned_pdb, '} type {pdb} first 0 last -1 step 1 waitfor 1\n']))
        f_out.write('animate style Once\n')

        # set protein representations
        f_out.write('mol addrep 0\n')
        f_out.write('mol modselect 1 0 protein\n')
        f_out.write('mol modstyle 1 0 NewCartoon 0.300000 10.000000 4.100000 0\n')
        f_out.write('mol modcolor 1 0 ColorID 0\n')

        # set ligand representations
        f_out.write('mol addrep 0\n')
        f_out.write(''.join(['mol modselect 2 0 resname ', ligand, '\n']))
        f_out.write('mol modstyle 2 0 CPK 1.000000 0.300000 12.000000 12.000000\n')
        f_out.write('mol modcolor 2 0 ColorID 4\n')

        # set protein representations
        f_out.write('mol addrep 1\n')
        f_out.write('mol modselect 1 1 protein\n')
        f_out.write('mol modstyle 1 1 NewCartoon 0.300000 10.000000 4.100000 0\n')
        f_out.write('mol modcolor 1 1 ColorID 1\n')

        # set ligand representations
        f_out.write('mol addrep 1\n')
        f_out.write(''.join(['mol modselect 2 1 resname ', ligand, '\n']))
        f_out.write('mol modstyle 2 1 CPK 1.000000 0.300000 12.000000 12.000000\n')
        f_out.write('mol modcolor 2 1 ColorID 7\n')

        # hide original representations
        f_out.write('mol showrep 0 0 0\n')
        f_out.write('mol showrep 1 0 0\n')



base_dir = '/home/mrama/Desktop/MD/GAPDH/2_Dock_March/1_APO/1_G3P/dock_prep/'
output_dir = '/home/mrama/Desktop/MD/GAPDH/2_Dock_March/1_APO/1_G3P/'

ref_number = str(649)
ref = ''.join([base_dir, 'GAPDH_WT_APO-', ref_number, '_receptor_G3P_docked.pdb'])

ligand = 'G3H'
frames = [308, 320, 153]

for frame_number in frames:
    frame_number = str(frame_number)
    frame = ''.join([base_dir, 'GAPDH_WT_APO-', frame_number, '_receptor_G3P_docked.pdb'])

    aligned_pdb = ''.join([output_dir, 'frame_', frame_number, '_aligned_to_ref_', ref_number, '.pdb'])
    align_enzymes(ref, frame, aligned_pdb)

    output_file = ''.join([output_dir, 'rmsd_frame_', frame_number, '_aligned_to_', ref_number])
    calculate_traj_rmsd_fit(ref, frame, output_file, 'backbone')

    file_out = ''.join([output_dir, 'viz_frame_', frame_number, '_aligned_to_', ref_number, '.tcl'])
    generate_vmd_viz(file_out, ref, aligned_pdb, ligand)
