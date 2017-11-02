from matplotlib import *
use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import re


def harvest_mmpbsa_perframe(base_dir, n_frames, file_out):
    """
    Gets the DeltaG from each frame, writes these in a csv file and plots it.

    :param base_dir: path to directory where each 'frame folder' is
    :param n_frames: how many frames were used
    :param file_out: file name where the results should be written/plotted
    :return: None
    """

    deltaG_all = {}

    for i in range(1, n_frames+1):
        print ''.join([base_dir, '/frame_', str(i), '/FINAL_RESULTS_MMPBSA.dat'])
        with open(''.join([base_dir, '/frame_', str(i), '/FINAL_RESULTS_MMPBSA.dat']), 'r') as f_in:
            content = f_in.read()
            deltaG_i = re.findall('DELTA\s+TOTAL\s+(\S+)\s+', content)[0]
            deltaG_all[i] = [float(deltaG_i)]

    deltaG_df = pd.DataFrame.from_dict(deltaG_all, orient='index')
    deltaG_df.index.name = 'frame'
    deltaG_df.columns = ['deltaG']

    deltaG_df.plot()
    plt.xlabel('frame number')
    plt.ylabel('Delta G')
    plt.savefig(''.join([base_dir, '/', file_out, '.pdf']))
    plt.savefig(''.join([base_dir, '/', file_out, '.png']))
    plt.close()

    dG_avg = deltaG_df.mean(axis=0)
    dG_std = deltaG_df.std(axis=0)

    deltaG_df.loc[len(deltaG_df.index) + 1] = dG_avg.values
    deltaG_df.loc[len(deltaG_df.index) + 2] = dG_std.values

    deltaG_df.to_csv(''.join([base_dir, '/', file_out, '.dat']), sep='\t')


# TODO: implement error propagation
def average_dG(base_dir, clusters, file_out):
    """
    For each enzyme + ligand pair, goes to each cluster and averages the deltaG values over all clusters.

    :param base_dir: path to where the 'cluster folders' are
    :param clusters: a list of the cluster folders, e.g. ['cl001', 'cl111', 'cl010']
    :param file_out: file name where the results should be written
    :return: None
    """

    deltaG_all = {}

    for cl in clusters:
        with open('/'.join([base_dir, cl, 'FINAL_RESULTS_MMPBSA.dat']), 'r') as f_in:
            content = f_in.read()
            try:
                deltaG_i = re.findall('DELTA\s*TOTAL\s+(\S+)\s+(\S+)\s+(\S+)', content)[0]
            except IndexError:
                raise IndexError('The line on DELTA TOTAL wasn\'t found.')
            deltaG_all[cl] = [float(deltaG_i[0]), float(deltaG_i[1]), float(deltaG_i[2])]

    deltaG_df = pd.DataFrame.from_dict(deltaG_all, orient='index')
    deltaG_df.index.name = 'frame'
    deltaG_df.columns = ['deltaG', 'std', 'sem']

    # append mean and std values
    dG_avg = deltaG_df.mean(axis=0)
    dG_std = deltaG_df.std(axis=0)

    deltaG_df.loc[len(clusters)] = dG_avg.values
    deltaG_df.loc[len(clusters) + 1] = dG_std.values
    deltaG_df = deltaG_df.rename(index={len(clusters): 'mean', len(clusters) + 1: 'std'})

    deltaG_df.to_csv(''.join([base_dir, '/', file_out, '.dat']), sep='\t')


# TODO: implement error propagation
def get_ddG(base_dir, substrate_dir, product_dir, clusters, file_out):
    """
    For each enzyme + (substrate, product) calculates the difference in binding energies between binding the
    product and binding the substrate, resulting in a DeltaDeltaG.

    :param base_dir: path to where the enzyme form is, e.g. '/home/user/TALB/1_APO'
    :param substrate_dir: folder where the substrate related results are, e.g. '1_S7P'
    :param product_dir: folder where the product related results are, e.g. '2_F6P'
    :param file_out: file name where the results should be written
    :return: None
    """

    substrate_df = pd.read_csv(''.join([base_dir, substrate_dir, '/dG_sumup.dat']), sep='\t')
    deltaG_subs = substrate_df['deltaG'][len(clusters)]

    product_df = pd.read_csv(''.join([base_dir, product_dir, '/dG_sumup.dat']), sep='\t')
    deltaG_prod = product_df['deltaG'][len(clusters)]

    ddG = round(deltaG_prod - deltaG_subs, ndigits=3)

    with open(''.join([base_dir, '/', file_out, '.dat']), 'w') as f_out:
        f_out.write(str(ddG))


def scp_mmpbsa_res(server, from_folder, to_folder):
    """
        Copy out, crd, prmtop files from server to local machine

    :param server: server name, e.g. ssb3.ucsd.edu
    :param from_folder: where the files are located
    :param to_folder: where the files should be copied to
    """

    from paramiko import SSHClient
    from scp import SCPClient

    ssh = SSHClient()
    ssh.load_system_host_keys()
    ssh.connect(server, username='marta')

    with SCPClient(ssh.get_transport(), sanitize=lambda x: x) as scp:
        scp.get(''.join([from_folder, '/*.dat']), to_folder, preserve_times=True, recursive=True)
        scp.get(''.join([from_folder, '/*.pdf']), to_folder, preserve_times=True, recursive=True)
        scp.get(''.join([from_folder, '/*.png']), to_folder, preserve_times=True, recursive=True)


def run_harvest_per_frame(enzyme_dir, clusters, file_out, n_frames_list, per_frame_folder):

    base_dir = ''.join(['/home/marta/', enzyme_dir])

    for i, cl in enumerate(clusters):
        cl_dir = '/'.join([base_dir, cl, per_frame_folder])
        cl_file_out = '_'.join([file_out, cl])
        harvest_mmpbsa_perframe(cl_dir, n_frames_list[i], cl_file_out)


def copy_to_local_machine(desktop_dir, enzyme_dir, clusters, server, per_frame_folder):

    base_dir = ''.join(['/home/marta/', enzyme_dir])

    for i, cl in enumerate(clusters):
        cl_dir = '/'.join([base_dir, cl, per_frame_folder])

        to_folder = ''.join([desktop_dir, enzyme_dir])
        print to_folder
        if not os.path.isdir(to_folder):
            os.makedirs(to_folder)
        scp_mmpbsa_res(server, cl_dir, to_folder)


# to delete
def run_harvest():

    base_dir = '/home/marta/ENO_1E9I/4_MMPBSA2016/1_APO/2_PEP'
    clusters = ['cl110_2', 'cl100_2', 'cl000_2']
    file_out = 'ENO_1E9I_WT_APO_PEP_2_dG_sumup'
    average_dG(base_dir, clusters,  file_out)


if __name__ == '__main__':

    enzyme_dir = 'ENO_AB/4_MMPBSA2016/1_MG/2_PEP'
    clusters = ['cl000', 'cl302', 'cl312', 'cl402']
    file_out = 'ENO_AB_WT_APO_PEP'
    n_frames_list = [500, 140, 140, 140]
    per_frame_folder = 'per_frame'

    #run_harvest_per_frame(enzyme_dir, clusters, file_out, n_frames_list, per_frame_folder)

    server = '137.110.115.63'
    desktop_dir = '/home/mrama/Desktop/MD/'
    copy_to_local_machine(desktop_dir, enzyme_dir, clusters, server, per_frame_folder)