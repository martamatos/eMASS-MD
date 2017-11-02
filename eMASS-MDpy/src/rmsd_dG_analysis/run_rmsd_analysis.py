import numpy as np
import pandas as pd
from collections import OrderedDict
import matplotlib.pyplot as plt

from src.kinetics_integration.plots_definitions import plot_timecourses_defs


def gather_rmsd_plots(root_dir, form_list, ligand_list, cluster_list, n_frames_per_form_per_ligand, file_out):
    """

    Plot all rmsd values for each MD simulation in the same figure.

    :param root_dir: path to folder with all MD simulations.
    :param form_list: list with enzyme forms considered, e.g., ['APO', 'HALO_S7P_remS7P']
    :param ligand_list: list with ligands considered per enzyme form  ['2PG', 'PEP']
    :param cluster_list: list of cluster per enzyme form per ligand
    :param n_frames_per_form_per_ligand: how many frames per enzyme form per ligand in the simulation
    :param file_out: path+name for output plot
    :return: None
    """

    form_folder_list = ['_'.join([str(i+1), form]) for i, form in enumerate(form_list)]
    ligand_folder_list = ['_'.join([str(i+1), ligand]) for i, ligand in enumerate(ligand_list)]

    n_rows = len(form_list)*n_frames_per_form_per_ligand
    fig, ax = plt.subplots(nrows=n_rows, ncols=len(ligand_list)*2, figsize=[25, n_rows*3], sharex=False, sharey=True)
    min_rmsd = 10**10
    max_rmsd = -10**10

    for ligand_i in range(len(ligand_list)):
        for form_i in range(len(form_list)):
            for cluster_i, cluster in enumerate(cluster_list[form_i][ligand_i]):
                for group_i, group in enumerate([ligand_list[ligand_i], 'backbone_100']):

                    if group == 'S7P_linear':
                        group = 'S7P'

                    base_dir = ''.join([root_dir, form_folder_list[form_i], '/', ligand_folder_list[ligand_i], '/', cluster, '/'])
                    file_name = ''.join(['rmsd_', group, '_', cluster])
                    file_in = ''.join([base_dir, file_name, '.dat'])

                    data_df = pd.read_csv(file_in, sep=' ', index_col=0, header=None)
                    rmsd_values = data_df.ix[:, 2].values

                    if max(rmsd_values) > max_rmsd:
                        max_rmsd = max(rmsd_values)
                    if min(rmsd_values) < min_rmsd:
                        min_rmsd = min(rmsd_values)

                    ax[form_i*n_frames_per_form_per_ligand + cluster_i][ligand_i*2 + group_i].plot(rmsd_values)
                    ax[form_i*n_frames_per_form_per_ligand + cluster_i][ligand_i*2 + group_i].set_title(', '.join([ligand_list[ligand_i], form_list[form_i],  cluster, group]))

    for i in range(len(ax)):
        for j in range(len(ax[0])):
            ax[i][j].set_ylim(min_rmsd, max_rmsd)
            ax[len(ax)-1][j].set_xlabel('frame number')
        ax[i][0].set_ylabel('rmsd (A)')

    plt.tight_layout()
    plt.savefig(''.join([file_out, '.pdf']))
    plt.close()


def plot_rmsd_values_scatter(data_df, file_out, substrate, product):
    """

    Plots rmsd values average+-std in a scatter plot.

    :param data_df: pandas dataframe with all rmsd values
    :param file_out: path+name for output plot
    :param substrate: substrate name
    :param product: product name
    :return: None
    """

    plot_timecourses_defs()

    rmsd_substrate_avg_list = data_df[(data_df['ligand'] == substrate) & (data_df['group'] == substrate)]['mean'].values
    rmsd_substrate_std_list = data_df[(data_df['ligand'] == substrate) & (data_df['group'] == substrate)]['std'].values
    rmsd_backbone_substrate_avg_list = data_df[(data_df['ligand'] == substrate) & (data_df['group'] == 'backbone')]['mean'].values
    rmsd_backbone_substrate_std_list = data_df[(data_df['ligand'] == substrate) & (data_df['group'] == 'backbone')]['std'].values

    rmsd_product_avg_list = data_df[(data_df['ligand'] == product) & (data_df['group'] == product)]['mean'].values
    rmsd_product_std_list = data_df[(data_df['ligand'] == product) & (data_df['group'] == product)]['std'].values
    rmsd_backbone_product_avg_list = data_df[(data_df['ligand'] == product) & (data_df['group'] == 'backbone')]['mean'].values
    rmsd_backbone_product_std_list = data_df[(data_df['ligand'] == product) & (data_df['group'] == 'backbone')]['std'].values

    xtick_labels_substrate = list(data_df[(data_df['ligand'] == substrate) & (data_df['group'] == substrate)].index.values)
    xtick_labels_substrate.insert(0, '')
    xtick_labels_product = list(data_df[(data_df['ligand'] == product) & (data_df['group'] == product)].index.values)
    xtick_labels_product.insert(0, '')


    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8,6), sharex=True, sharey=False)

    ax[0, 0].errorbar(x=range(len(rmsd_backbone_substrate_avg_list)), y=rmsd_backbone_substrate_avg_list, yerr=rmsd_backbone_substrate_std_list, fmt='o',  color='#377eb8', markersize=8)
    ax[1, 0].errorbar(x=range(len(rmsd_substrate_avg_list)), y=rmsd_substrate_avg_list, yerr=rmsd_substrate_std_list, fmt='o',  color='#377eb8', markersize=8)

    ax[0, 1].errorbar(x=range(len(rmsd_backbone_product_avg_list)), y=rmsd_backbone_product_avg_list, yerr=rmsd_backbone_product_std_list, fmt='o',  color='#377eb8', markersize=8)
    ax[1, 1].errorbar(x=range(len(rmsd_product_avg_list)), y=rmsd_product_avg_list, yerr=rmsd_product_std_list, fmt='o',  color='#377eb8', markersize=8)

    if substrate == 'S7P_linear':
        substrate = 'S7P'

    ax[0, 0].set_title(substrate + ' backbone')
    ax[0, 1].set_title(product + ' backbone')
    ax[1, 0].set_title(substrate + ' ligand')
    ax[1, 1].set_title(product + ' ligand')

    for i in range(len(ax)):
        for j in range(len(ax[0])):
            ax[i, j].set_xlim(-1, len(rmsd_substrate_avg_list))
            ax[i, j].set_ylim(0, 2.6)
            ax[i, j].yaxis.grid(True)
            ax[1, j].set_xlabel('measurement')

        ax[i, 0].set_ylabel('rmsd ($\AA$)')

    ax[0, 1].set_xlim(-1, len(rmsd_product_avg_list))


    plt.sca(ax[1, 0])  # set axis current axis active
    plt.xticks(range(-1, len(rmsd_substrate_avg_list)), xtick_labels_substrate, rotation=90)
    plt.sca(ax[1, 1])  # set axis current axis active
    plt.xticks(range(-1, len(rmsd_product_avg_list)), xtick_labels_product, rotation=90)

    plt.tight_layout()
    plt.savefig(''.join([file_out, '.pdf']))
    plt.close()



def sumup_rmsd(form_list, ligand_list, cluster_list, root_dir, file_out):
    """
    Sums up all MD trajectories rmsd results in a single file.

    :param form_list: list with enzyme forms considered, e.g., ['APO', 'HALO_S7P_remS7P']
    :param ligand_list: list with ligands considered per enzyme form  ['2PG', 'PEP']
    :param cluster_list: list of cluster per enzyme form per ligand
    :param root_dir: path to folder with all MD simulations
    :param file_out: file where rmsd results will be stored
    :return: None
    """

    form_folder_list = ['_'.join([str(i+1), form]) for i, form in enumerate(form_list)]
    ligand_folder_list = ['_'.join([str(i+1), ligand]) for i, ligand in enumerate(ligand_list)]

    rmsd_dic = OrderedDict()
    for form_i in range(len(form_list)):
        for ligand_i in range(len(ligand_list)):
            for cluster_i, cluster in enumerate(cluster_list[form_i][ligand_i]):
                for group_i, group in enumerate([ligand_list[ligand_i], 'backbone_100']):

                    if group == 'S7P_linear':
                        group = 'S7P'

                    base_dir = ''.join([root_dir, form_folder_list[form_i], '/', ligand_folder_list[ligand_i], '/', cluster, '/'])
                    file_name = ''.join(['rmsd_', group, '_', cluster])
                    file_in = ''.join([base_dir, file_name, '.dat'])

                    if group == 'backbone_100':
                        group = 'backbone'
                    if group == 'S7P':
                        group = 'S7P_linear'

                    data_df = pd.read_csv(file_in, sep=' ', index_col=0, header=None)
                    rmsd_values = data_df.ix[:, 2].values

                    rmsd_dic['_'.join([form_list[form_i], ligand_list[ligand_i], group, cluster])] = [ligand_list[ligand_i], group, np.mean(rmsd_values), np.var(rmsd_values), np.std(rmsd_values)]

    data_df = pd.DataFrame.from_dict(rmsd_dic, orient='index')
    data_df.columns = ['ligand', 'group', 'mean', 'var', 'std']

    data_df['std %'] = abs(data_df['std'] / data_df['mean'] * 100)

    data_df.to_csv(''.join([file_out, '.csv']), sep='\t')
    plot_rmsd_values_scatter(data_df, file_out, ligand_list[0], ligand_list[1])


def analyze_ENO(base_dir):

    root_dir = ''.join([base_dir, '/ENO_1E9I/3_MD_post_dock2016/'])
    file_out = ''.join([base_dir, '/ENO_1E9I/3_MD_post_dock2016/ENO_rmsd_sumup'])
    form_list = ['APO']
    ligand_list = ['2PG', 'PEP']
    cluster_list = [[['cl000', 'cl000_2', 'cl000_3', 'cl000_4', 'cl100', 'cl100_2', 'cl100_3', 'cl100_4', 'cl200', 'cl200_2', 'cl200_3', 'cl200_4'],
                    ['cl000', 'cl000_2', 'cl000_3', 'cl000_4', 'cl100', 'cl100_2', 'cl100_3', 'cl100_4', 'cl110', 'cl110_2', 'cl110_3', 'cl110_4']]]

    sumup_rmsd(form_list, ligand_list, cluster_list, root_dir, file_out)

    file_out = ''.join([base_dir, '/ENO_1E9I/3_MD_post_dock2016/ENO_all_rmsd'])
    n_frames_per_form_per_ligand = 12
    gather_rmsd_plots(root_dir, form_list, ligand_list, cluster_list, n_frames_per_form_per_ligand, file_out)


def analyze_ENO_AB(base_dir):

    root_dir = ''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/'])
    file_out = ''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/ENO_rmsd_sumup'])
    form_list = ['MG']
    ligand_list = ['2PG', 'PEP']
    cluster_list = [[['cl000', 'cl102', 'cl202', 'cl203'],
                    ['cl000', 'cl302', 'cl312', 'cl402']]]

    sumup_rmsd(form_list, ligand_list, cluster_list, root_dir, file_out)

    file_out = ''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/ENO_all_rmsd'])
    n_frames_per_form_per_ligand = 4
    gather_rmsd_plots(root_dir, form_list, ligand_list, cluster_list, n_frames_per_form_per_ligand, file_out)

    # this file is a modified version of enzyme_name_rmsd_sumup where the first column was changed, everything else is the same
    data_df = pd.read_csv(''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/ENO_rmsd_sumup_paper.csv']), sep=',', index_col=0)
    substrate = '2PG'
    product = 'PEP'
    file_out = ''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/ENO_rmsd_paper'])
    plot_rmsd_values_scatter(data_df, file_out, substrate, product)


def analyze_GAPD(base_dir):

    root_dir = ''.join([base_dir, '/GAPDH/3_MD_post_dock2016/'])
    file_out = ''.join([base_dir, '/GAPDH/3_MD_post_dock2016/GAPD_rmsd_sumup'])
    form_list = ['NAD', 'NAD_remG3P']
    ligand_list = ['G3P', 'DPG']
    cluster_list = [[['cl000', 'cl000_2', 'cl001', 'cl001_2', 'cl012', 'cl012_2'],
                     ['cl100', 'cl100_2', 'cl101', 'cl101_2', 'cl220', 'cl220_2']],
                    [['cl010', 'cl010_2', 'cl100', 'cl100_2', 'cl200', 'cl200_2'],
                     ['cl100', 'cl100_2', 'cl120', 'cl120_2', 'cl210', 'cl210_2']]]


    sumup_rmsd(form_list, ligand_list, cluster_list, root_dir, file_out)

    file_out = ''.join([base_dir, '/GAPDH/3_MD_post_dock2016/GAPD_all_rmsd'])
    n_frames_per_form_per_ligand = 6
    gather_rmsd_plots(root_dir, form_list, ligand_list, cluster_list, n_frames_per_form_per_ligand, file_out)

    # this file is a modified version of enzyme_name_rmsd_sumup where the first column was changed, everything else is the same
    data_df = pd.read_csv(''.join([base_dir, '/GAPDH/3_MD_post_dock2016/GAPD_rmsd_sumup_paper.csv']), sep=',', index_col=0)
    substrate = 'G3P'
    product = 'DPG'
    file_out = ''.join([base_dir, '/GAPDH/3_MD_post_dock2016/GAPD_rmsd_paper'])
    plot_rmsd_values_scatter(data_df, file_out, substrate, product)


def analyze_TALB(base_dir):

    root_dir = ''.join([base_dir, '/TalB/3_MD_post_dock2016/'])
    file_out = ''.join([base_dir, '/TalB/3_MD_post_dock2016/TALB_rmsd_sumup'])
    form_list = ['APO', 'HALO_S7P_remS7P']
    ligand_list = ['S7P_linear', 'F6P']
    cluster_list = [[['cl101', 'cl101_2', 'cl102', 'cl102_2', 'cl212', 'cl212_2'],
                     ['cl202', 'cl202_2', 'cl212', 'cl212_2', 'cl321', 'cl321_2']],
                    [['cl101', 'cl101_2', 'cl302', 'cl302_2', 'cl322', 'cl322_2'],
                     ['cl101', 'cl101_2', 'cl201', 'cl201_2', 'cl311', 'cl311_2']]]
    sumup_rmsd(form_list, ligand_list, cluster_list, root_dir, file_out)

    file_out = ''.join([base_dir, '/TalB/3_MD_post_dock2016/TALB_all_rmsd'])
    n_frames_per_form_per_ligand = 6
    gather_rmsd_plots(root_dir, form_list, ligand_list, cluster_list, n_frames_per_form_per_ligand, file_out)

    # this file is a modified version of enzyme_name_rmsd_sumup where the first column was changed, everything else is the same
    data_df = pd.read_csv(''.join([base_dir, '/TalB/3_MD_post_dock2016/TALB_ddG_sumup_paper.csv']), sep=',', index_col=0)
    substrate = 'S7P_linear'
    product = 'F6P'
    file_out = ''.join([base_dir, '/TalB/3_MD_post_dock2016/TALB_rmsd_paper'])
    plot_rmsd_values_scatter(data_df, file_out, substrate, product)


def main():
    base_dir = '/home/mrama/Dropbox/PhD_stuff/Projects/MD/eMASS-MD/MD_data'

    analyze_ENO(base_dir)
    analyze_ENO_AB(base_dir)
    analyze_GAPD(base_dir)
    analyze_TALB(base_dir)


if __name__ == '__main__':
    main()