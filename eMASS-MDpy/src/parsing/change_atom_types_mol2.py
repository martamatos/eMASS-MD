import re


def change_atom_types_mol2(filename, ligand, atoms):
    """
    Given a mol2 file, changes the atom types.

    :param filename: path+name for mol2 file
    :param ligand: ligand name
    :param atoms: list of atoms if these are to be substituted
    :return: None
    """

    try:
        f_in = open(filename + ".mol2", "r")
    except:
        print "oops", filename
        return

    # file where the modified mol2 content will be written into
    f_out = open(filename + "2.mol2", "w+")

    i = 0
    line = f_in.readline()
    while line:

        if re.findall(" (\w+)\s*1[ ]{2}" + ligand, line):
            # so that the column after the atom type starts in the right place
            sp = 8 - len(atoms[i])
            # if there is an atom type + charge in the current line, substitute the atom
            #  type by the on in new_atom_list and leave the charge with only 4 decimal
            #  places
            res = re.sub(" (\w+)\s*1[ ]{2}" + ligand, " " + atoms[i] + " "*sp + "1  " + ligand, line)

            i += 1
        else:
            i = 0
            res = line
        f_out.write(res)
        line = f_in.readline()

    f_in.close()
    f_out.close()


new_atom_list = ["H", "C.3", "C.3", "C.3", "O.3", "P.3", "O.co2", "C.3", "C.3", "O.3",
                 "O.co2", "O.co2", "O.3", "H", "H", "O.3", "H", "H", "H", "H", "H", "H",
                 "C.2", "C.3", "H", "C.3", "O.2", "O.3", "H", "H", "H"]

#n_frames = 803
#ligand = "S7P"
#for i in range(1, n_frames + 1):
#    filename = "/home/marta/Dock/enzymes/TALB/1_APO/1_S7P/dock_prep/TALB_WT_APO-%.3d_S7P_flexible_scored" % (i)
#    #filename = "../tests/test_files/change_atom_types_mol2/TALB_WT_APO-%.3d_S7P_flexible_scored" % (i)
#    change_atom_types_mol2(filename, ligand, new_atom_list)