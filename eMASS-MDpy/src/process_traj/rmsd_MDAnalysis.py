import matplotlib.pyplot as plt
import pandas as pd
from MDAnalysis import Universe
from MDAnalysis.analysis.rms import RMSD


def calculate_traj_rmsd_fit(traj, prmtop, file_out, selection):
    """
        Calculates the RMSD of the protein's backbone in a given MD trajectory relative to its first frame.
        The trajectory is aligned first.
        Saves a file with the rmsd for each frame, as well as a plot.

    :param traj: path to the MD trajectory is.
    :param prmtop: path to the prmtop file.
    :param file_out: path to the file with output rmsd and plot WITHOUT THE FILE EXTENSION
    :param selection: for which part of the system will the RMSD be calculated, e.g. backbone, resname G6P
    :return: None
    """

    ref = Universe(prmtop, traj)
    ref.trajectory[0]  # go to first frame
    trj = Universe(prmtop, traj)

    R = RMSD(trj, ref, select=selection, filename='.'.join([file_out, 'dat']))
    R.run()
    R.save()

    rmsd = R.rmsd.T
    plt.plot(rmsd[0], rmsd[2], 'k-', linewidth=2)

    plt.xlabel('frame number')
    plt.ylabel(r'RMSD ($\AA$)')
    plt.savefig('.'.join([file_out, 'png']))
    plt.close()


def average_rmsd_from_file(file_rmsd, file_out):
    """
        Given a file with the rmsd, as writen by calculate_traj_rmsd_fit, averages the protein rmsd and
        calculates the std, and writes the result in rmsd_sumup.dat

    :param file_rmsd: file with the rmsd, where the rmsd is in the second column and columns are separated by one space.
    :param file_out: path to the file with output rmsd and plot WITHOUT THE FILE EXTENSION
    :return: None
    """

    df = pd.read_csv(file_rmsd, sep=' ', header=None)

    with open(file_out, 'w') as f_out:
        f_out.write(''.join(['mean:\t', str(df[2].mean()), '\n']))
        f_out.write(''.join(['std:\t', str(df[2].std()), '\n']))
