import re


def _get_atom_order_off(off_file, n_atoms):

    with open('.'.join([off_file, 'off']), 'r') as f_off:
        content = f_off.read()
    res = re.findall(' \"(\w+)\" \"\w+\"', content)

    return res[:n_atoms]


def _get_atom_order_pdb(pdb_file, ligand):
    map_atom_to_line = {}

    with open('.'.join([pdb_file, 'pdb']), 'r') as f_pdb:
        line = f_pdb.readline()
        while line:
            #res = re.findall(''.join(['ATOM\s+\d+\s*(\w+)\s+', ligand]), line)
            res = re.findall(''.join(['HETATM\s+\d+\s*(\w+)\s+', ligand]), line)
            if res:
                map_atom_to_line[res[0]] = line
            line = f_pdb.readline()

    return map_atom_to_line


def _write_ligand(f_pdb_out, atom_order_off, map_atom_to_line):
    for atom in atom_order_off:
        f_pdb_out.write(map_atom_to_line[atom])


# TODO generalize for multiple ligands
def reorder_ligand_atoms_from_off(pdb_file, off_file, ligand, n_atoms):
    """
    Given a pdb file and respective off file, it reorders the ligand atoms in the pdb file according to the respective
    off file.

    :param pdb_file: path+name for the pdb file
    :param off_file: path+name for the off file
    :param ligand: ligand name
    :param n_atoms: number of atoms in the ligand
    :return: None
    """

    map_atom_order = _get_atom_order_off(off_file, n_atoms)
    map_atom_to_line = _get_atom_order_pdb(pdb_file, ligand)

    with open('.'.join([pdb_file, 'pdb']), 'r') as f_pdb, open(''.join([pdb_file, '2.pdb']), 'w+') as f_pdb_out:
        flag = 0
        line = f_pdb.readline()
        while line:

            if re.match(''.join(['HETATM\s+\d+\s+\w+\s+', ligand]), line):
            #if re.match(''.join(['ATOM\s+\d+\s+\w+\s+', ligand]), line):
                if not flag:
                    f_pdb_out.write('TER\n')
                    _write_ligand(f_pdb_out, map_atom_order, map_atom_to_line)
                    flag = 1
            else:
                f_pdb_out.write(line)

            line = f_pdb.readline()


def rename_ligand_atoms(pdb_file, ligand_list, map_atom_names=None):
    """
    Given a PDB file where the ligand is a 'HETATM' and not a 'ATOM', renames the ligand atoms.
    If no dictionary is given with the atom names mapping, the atom names will be updated as follows:
      - 'C' -> 'C1'
      - 'O1' -> 'O2'

    Otherwise they will be updated according to map_atom_names.

    :param pdb_file: file provided to used to rename the atom
    :param ligand_list: a list of ligands whose atoms need to be renamed
    :param map_atom_names: a dictionary mapping the ligand atom names in pdb_file to the new ones.
    :return: None
    """

    for ligand in ligand_list:

        with open('.'.join([pdb_file, 'pdb']), 'r') as f_pdb, open(''.join([pdb_file, '2.pdb']), 'w+') as f_pdb_out:
            line = f_pdb.readline()
            while line:
                white_spaces = 2  # to deal with spaces after the atom name

                ligand_match = re.findall(''.join(['HETATM\s+\d+\s+\w+\s+(', ligand, ')']), line)

                if ligand_match:
                    if map_atom_names:
                        res = re.findall(''.join(['HETATM\s+\d+\s+([A-Z]+[0-9]*)(\s+)', ligand]), line)[0]
                        atom = res[0]
                        white_spaces = len(res[1]) - (len(map_atom_names[atom]) - len(atom))
                        new_line = re.sub(''.join(['(HETATM\s+\d+\s+)([A-Z]+[0-9]*)(\s+', ligand, ')']), ''.join([r'\1', map_atom_names[atom], ' '*white_spaces, ligand]) , line)

                    else:
                        atom_parts = re.findall(''.join(['HETATM\s+\d+\s+([A-Z]+)([0-9]*)\s+', ligand]), line)[0]
                        if atom_parts[-1] == '':
                            atom_sub = ''.join([atom_parts[0], '1'])
                            white_spaces = 2
                        else:
                            new_atom_number = str(int(atom_parts[1]) + 1)
                            atom_sub = ''.join([atom_parts[0], new_atom_number])
                            if len(new_atom_number) == 2:
                                white_spaces = 1

                        new_line = re.sub(''.join(['(HETATM\s+\d+\s+)([A-Z]+[0-9]*)(\s+', ligand, ')']), ''.join([r'\1', atom_sub, ' '*white_spaces, ligand]) , line)
                    f_pdb_out.write(new_line)
                else:
                    f_pdb_out.write(line)

                line = f_pdb.readline()


