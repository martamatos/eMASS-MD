import re


def parse_vmd_rmsd_output(rmsd_file):
    """
    Given a trajectory rmsd file created through vmd, transform it into a proper TSV file.

    :param rmsd_file: file path+name for rmsd file.
    :return: None
    """

    with open(rmsd_file, 'r') as f_in, open(rmsd_file[:-4] + '_mod.dat', 'w') as f_out:
        line = f_in.readline()
        while line:
            res = re.sub('^\s*(\S+)\s+(\S+)', r'\1\t\2', line)
            f_out.write(res)
            line = f_in.readline()


#parse_vmd_rmsd_output('../tests/test_files/calculate_rmsd_MDanalysis/true_res_trajrmsd_dcd.dat')