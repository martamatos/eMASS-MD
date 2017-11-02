import re


def renumber_pdb_atoms(file_in, file_out):
    """
    Given a pdb file it renumbers the atoms to start in 1.

    :param file_in: path+name of input pdb file
    :param file_out: path+name of output pdb file
    :return: None
    """

    atom_i = 1

    with open(file_in, 'r') as f_in,  open(file_out, 'w') as f_out:
        line = f_in.readline()

        while line:
            sub = re.sub('(ATOM\s+)(\d+)(.*)\n', ''.join([r'\1 ', str(atom_i), r'\3\n']), line)
            f_out.write(sub)

            line = f_in.readline()
            atom_i += 1


file_in = '/home/mrama/Desktop/MD/TalB/0_Crystalline_structures/4s2c_A.pdb'
file_out = '/home/mrama/Desktop/MD/TalB/0_Crystalline_structures/4s2c_A_mod.pdb'
renumber_pdb_atoms(file_in, file_out)