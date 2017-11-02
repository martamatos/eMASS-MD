import re


def substitute_ligands(pdb_file, old_ligand, new_ligand, atom_map):
    """
    Given a pdb file, substitute the old_ligand by new_ligand.
    Seems to assume same number of atoms.

    :param pdb_file: path+name for pdb_file
    :param old_ligand: old ligand name
    :param new_ligand: new ligand name
    :param atom_map: dictionary mapping old ligand atoms to new ligand atoms
    :return: None
    """

    pdb_file_out = ''.join([pdb_file, '2.pdb'])

    with open(''.join([pdb_file, '.pdb']), 'r') as f_in, open(pdb_file_out, 'w+') as f_out:

        line = f_in.readline()
        while line:

            if re.search(old_ligand, line):

                while re.search(old_ligand, line):
                    print line
                    old_atom = re.findall('HETATM\s+\d+\s+([A-Z0-9]+[0-9]*)', line)[0]
                    new_atom = atom_map[old_atom] if old_atom in atom_map.keys() else old_atom

                    res = re.sub(''.join(['(HETATM\s+\d+\s+)([A-Z0-9]+[0-9]*)(\s+)(', old_ligand, ')(.*\n)']),
                                 ''.join([r'\g<1>', new_atom, r'\g<3>', new_ligand, r'\g<5>']),
                                 line)
                    print res
                    print '----------------'
                    f_out.write(res)
                    line = f_in.readline()
                f_out.write(line)
            else:
                f_out.write(line)

            line = f_in.readline()


if __name__ == '__main__':

    #pdb_file = '/home/mrama/Desktop/MD/GAPDH/1_MD/3_HALO/GAPDH_WT_HALO_G3P_lasframe_noWATNa'
    #old_ligand = 'G3H'
    #new_ligand = 'DPG'
    #atom_map = {'P1': 'P2', 'O3': 'O8', 'O4': 'O9', 'O5': 'O10', 'C1': 'C3', 'C2': 'C2', 'C3': 'C1', 'H1': 'H3',
    #            'H2': 'H4', 'O2': 'O6', 'H3': 'H1', 'H5': 'H2', 'O6': 'O5'}

    #pdb_file = '/home/mrama/Desktop/MD/GAPDH/1_MD/3_HALO/GAPDH_WT_HALO_G3P_lasframe_noWATNa_G3P'
    #old_ligand = 'NAD'
    #new_ligand = 'NAD'
    #atom_map = {'H86': 'H73', 'H87': 'H67', 'H92': 'H66', 'H81': 'H71', 'H89': 'H70', 'H90': 'H68',
    #            'H77': 'H62', 'H79': 'H64', 'H76': 'H61', 'H78': 'H63', 'H75': 'H60', 'H80': 'H65',
    #            'H73': 'H58', 'H74': 'H59',
    #            'H52': 'H48', 'H51': 'H47',
    #            'H53': 'H49', 'H55': 'H51', 'H57': 'H53', 'H54': 'H50', 'H56': 'H52', 'H58': 'H54',
    #            'H60': 'H55', 'H67': 'H57', 'H65': 'H56', 'H64': 'H72'}


    pdb_file = '/home/mrama/Desktop/MD/G6PD/1_MD/2_NADP_single/G6PD_WT_NADPH_lasframe_6PG_edit'
    old_ligand = 'G6P'
    new_ligand = '6PG'
    atom_map = {'O6': 'O4', 'C6': 'C5', 'C3': 'C2', 'H6': 'H5', 'H7': 'H6',
                'C1': 'C1', 'C2': 'C3', 'C4': 'C4', 'C5': 'C6', 'O1': 'O1',
                'H1': 'H1', 'O2': 'O2', 'H8': 'H7',
                'O3': 'O3', 'H9': 'H8', 'H2': 'H3',
                'O4': 'O5', 'H10': 'H9', 'H4': 'H4',
                'O5': 'O6', 'H3': 'H2'}

    substitute_ligands(pdb_file, old_ligand, new_ligand, atom_map)



