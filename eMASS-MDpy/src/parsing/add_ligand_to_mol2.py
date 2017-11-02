import re


def add_info_to_mol2_file(filename, ligand, atoms=None):
    """
    Given a mol2 file where the ligand name is missing, it ads the ligand name.
    If atoms is provided as a list of new atoms, the ligand atoms will also be substituted by the new ones.

    :param filename: path+name for mol2 file
    :param ligand: ligand name
    :param atoms: list of atoms if these are to be substituted (optional)
    :return: None
    """

    try:
        f_in = open(filename + ".mol2", "r")
    except:
        return

    # file where the modified mol2 content will be written into
    f_out = open(filename + "2.mol2", "w+")

    i = 0
    line = f_in.readline()
    while line:

        # find the atom type (group 1) and the atom charge (group 2)
        charge = re.findall(" (\w+.?\w*)\s*1[ ]{2}([-]?[0-9].[0-9]{6})[ ]{4}0.0000", line)

        if charge:

            if atoms:
                # so that the column after the atom type starts in the right place
                sp = 8 - len(atoms[i])
                # get the charge, convert to float and back to string,
                #  then it has only 4 decimal places, not the best way to do it obviously, but works
                charge = str(float(charge[0][1]))

                # if there is an atom type + charge in the current line, substitute the atom
                #  type by the on in new_atom_list and leave the charge with only 4 decimal
                #  places

                res = re.sub(" (\w+.?\w*)\s*1[ ]{2}([-]?[0-9].[0-9]{6})[ ]{4}0.0000",
                             " " + atoms[i] + " "*sp + "1  " + ligand + r"         " + charge,
                             line)
            else:
                sp = 8 - len(charge[0][0])
                charge = str(float(charge[0][1]))
                res = re.sub(" (\w+.?\w*)\s*1[ ]{2}([-]?[0-9].[0-9]{6})[ ]{4}0.0000",
                             r' \1' + " "*sp + "1  " + ligand + r"         " + charge,
                             line)
            i += 1
        else:
            i = 0
            res = line
        f_out.write(res)
        line = f_in.readline()

    f_in.close()
    f_out.close()


# new_atom_list = ["H", "C.3", "C.3", "C.3", "O.3", "P.3", "O.co2", "C.3", "C.3", "O.3",
#                  "O.co2", "O.co2", "O.3", "H", "H", "O.3", "H", "H", "H", "H", "H", "H",
#                  "C.2", "C.3", "H", "C.3", "O.2", "O.3", "H", "H", "H"]
#
n_frames = 803
ligand = "S7P"
for i in range(1, n_frames + 1):
    filename = "/home/marta/Dock/enzymes/TALB/1_APO/1_S7P/dock_prep/TALB_WT_APO-%.3d_S7P_flexible_scored" % (i)
    #filename = "../tests/test_files/add_ligand_to_mol2/TALB_WT_APO-%.3d_S7P_flexible_scored" % (i)
    add_info_to_mol2_file(filename, ligand)