from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import shapiro
from scipy.stats import spearmanr

from src.kinetics_integration.plots_definitions import plot_timecourses_defs


def get_variance_of_the_average(variance_list):
    """

    Given a list with variance values, propagate variance for an average according to
        https://en.wikipedia.org/wiki/Propagation_of_uncertainty

    :param variance_list: a list with variance values.
    :return: the propagated variance
    """

    total_variance = sum(variance_list) / (len(variance_list)**2)

    return total_variance


def get_dKd_from_ddG(ddG):
    """

    Convert a ddG into a dKb

    :param ddG: ddG value in kcal/mol
    :return: None
    """

    R = 1.9872036 * 10**-3  # kcal K-1 mol-1
    T = 298.15  # K
    dKd = np.exp(ddG / (R*T))

    return dKd


def analyze_and_plot_dG_values_and_rmsd(data_df_rmsd, data_df_dG, file_out, substrate, product):
    """

    Plots dG values and rmsd values in different plots but same figure.
    (Need to run run_rmsd_analysis.py first)

    :param data_df_rmsd: pandas dataframe with rmsd data for each MD simulation
    :param data_df_dG: pandas daframe with dG data for each MD simulation
    :param file_out: path+name for plot file
    :param substrate: substrate bound to enzyme
    :param product: product bound to enzyme
    :return: None
    """

    rmsd_substrate_avg_list = data_df_rmsd[(data_df_rmsd['ligand'] == substrate) & (data_df_rmsd['group'] == substrate)]['mean'].values
    rmsd_substrate_std_list = data_df_rmsd[(data_df_rmsd['ligand'] == substrate) & (data_df_rmsd['group'] == substrate)]['std'].values
    rmsd_backbone_substrate_avg_list = data_df_rmsd[(data_df_rmsd['ligand'] == substrate) & (data_df_rmsd['group'] == 'backbone')]['mean'].values
    rmsd_backbone_substrate_std_list = data_df_rmsd[(data_df_rmsd['ligand'] == substrate) & (data_df_rmsd['group'] == 'backbone')]['std'].values

    rmsd_product_avg_list = data_df_rmsd[(data_df_rmsd['ligand'] == product) & (data_df_rmsd['group'] == product)]['mean'].values
    rmsd_product_std_list = data_df_rmsd[(data_df_rmsd['ligand'] == product) & (data_df_rmsd['group'] == product)]['std'].values
    rmsd_backbone_product_avg_list = data_df_rmsd[(data_df_rmsd['ligand'] == product) & (data_df_rmsd['group'] == 'backbone')]['mean'].values
    rmsd_backbone_product_std_list = data_df_rmsd[(data_df_rmsd['ligand'] == product) & (data_df_rmsd['group'] == 'backbone')]['std'].values

    xtick_labels_substrate = list(data_df_rmsd[(data_df_rmsd['ligand'] == substrate) & (data_df_rmsd['group'] == substrate)].index.values)
    xtick_labels_substrate.insert(0, '')
    xtick_labels_product = list(data_df_rmsd[(data_df_rmsd['ligand'] == product) & (data_df_rmsd['group'] == product)].index.values)
    xtick_labels_product.insert(0, '')

    fig = plt.figure(figsize=(10, 16))
    ax5 = fig.add_subplot(4, 2, 5)
    ax6 = fig.add_subplot(4, 2, 6)
    ax7 = fig.add_subplot(4, 2, 7)
    ax8 = fig.add_subplot(4, 2, 8)

    ax1 = fig.add_subplot(4, 2, 1, sharex=ax5)
    ax2 = fig.add_subplot(4, 2, 2, sharex=ax6)
    ax3 = fig.add_subplot(4, 2, 3, sharex=ax5)
    ax4 = fig.add_subplot(4, 2, 4, sharex=ax6)


    ax1.errorbar(x=range(len(rmsd_backbone_substrate_avg_list)), y=rmsd_backbone_substrate_avg_list, yerr=rmsd_backbone_substrate_std_list, fmt='o')
    ax3.errorbar(x=range(len(rmsd_substrate_avg_list)), y=rmsd_substrate_avg_list, yerr=rmsd_substrate_std_list, fmt='o')

    ax2.errorbar(x=range(len(rmsd_backbone_product_avg_list)), y=rmsd_backbone_product_avg_list, yerr=rmsd_backbone_product_std_list, fmt='o')
    ax4.errorbar(x=range(len(rmsd_product_avg_list)), y=rmsd_product_avg_list, yerr=rmsd_product_std_list, fmt='o')

    ax1.set_title('backbone')
    ax2.set_title('backbone')
    ax3.set_title('ligand')
    ax4.set_title('ligand')


    ax1.set_ylim(0, 3)
    ax2.set_ylim(0, 3)
    ax3.set_ylim(0, 3)
    ax4.set_ylim(0, 3)

    ax1.set_xlim(-1, len(rmsd_substrate_avg_list))
    ax2.set_xlim(-1, len(rmsd_substrate_avg_list))
    ax3.set_xlim(-1, len(rmsd_substrate_avg_list))
    ax4.set_xlim(-1, len(rmsd_substrate_avg_list))

    ax1.set_ylabel('rmsd (A)')
    ax3.set_ylabel('rmsd (A)')

    ax2.set_xlim(-1, len(rmsd_product_avg_list))

    plt.sca(ax3)  # set axis current axis active
    plt.xticks(range(-1, len(rmsd_substrate_avg_list)), (len(rmsd_product_avg_list)+1)*[''], rotation=90)
    plt.sca(ax4)  # set axis current axis active
    plt.xticks(range(-1, len(rmsd_product_avg_list)), (len(rmsd_product_avg_list)+1)*[''], rotation=90)

    plt.sca(ax1)  # set axis current axis active
    plt.xticks(range(-1, len(rmsd_substrate_avg_list)), (len(rmsd_product_avg_list)+1)*[''], rotation=90)
    plt.sca(ax2)  # set axis current axis active
    plt.xticks(range(-1, len(rmsd_product_avg_list)), (len(rmsd_product_avg_list)+1)*[''], rotation=90)

    #plt.tight_layout()
    #plt.savefig(''.join([file_out, '.pdf']))
    #plt.close()


    dG_substrate_list = data_df_dG[data_df_dG['ligand'] == substrate]['mean'][:-1].values
    dG_std_substrate_list = data_df_dG[data_df_dG['ligand'] == substrate]['std'][:-1].values

    dG_product_list = data_df_dG[data_df_dG['ligand'] == product]['mean'][:-1].values
    dG_std_product_list = data_df_dG[data_df_dG['ligand'] == product]['std'][:-1].values

    dG_mean_substrate = data_df_dG[data_df_dG['ligand'] == substrate]['mean'][-1]
    dG_mean_product = data_df_dG[data_df_dG['ligand'] == product]['mean'][-1]

    xtick_labels_substrate = list(data_df_dG[data_df_dG['ligand'] == substrate].index.values[:-1])
    xtick_labels_substrate.insert(0, '')
    xtick_labels_product = list(data_df_dG[data_df_dG['ligand'] == product].index.values[:-1])
    xtick_labels_product.insert(0, '')

    ymin = (np.min([np.min(dG_substrate_list-dG_std_substrate_list), np.min(dG_product_list-dG_std_product_list)])) - 5
    ymax = (np.max([np.max(dG_substrate_list+dG_std_substrate_list), np.max(dG_product_list+dG_std_product_list)])) + 5

    #fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10,7)) #, sharey=True)

    ax5.errorbar(x=range(len(dG_substrate_list)), y=dG_substrate_list, yerr=dG_std_substrate_list, fmt='o')
    ax5.axhline(y=dG_mean_substrate, color='black', linewidth=1.5)

    ax6.errorbar(x=range(len(dG_product_list)), y=dG_product_list, yerr=dG_std_product_list, fmt='o')
    ax6.axhline(y=dG_mean_product, color='black', linewidth=1.5)

    ax5.set_title(substrate)
    ax6.set_title(product)

    ax5.set_xlim(-1, len(dG_substrate_list))
    ax6.set_xlim(-1, len(dG_product_list))
    ax5.set_ylim(ymin, ymax)
    ax6.set_ylim(ymin, ymax)
    ax5.set_xlabel('measurement')
    ax6.set_xlabel('measurement')
    ax5.set_ylabel('dG (kcal/mol)')

    plt.sca(ax5)  # set axis current axis active
    xtick_labels_substrate = [ligand.replace('_linear_', '_') for ligand in xtick_labels_substrate]
    xtick_labels_substrate = [ligand.replace('HALO_S7P_', 'HALO_') for ligand in xtick_labels_substrate]
    xtick_labels_product = [ligand.replace('_linear_', '_') for ligand in xtick_labels_product]
    xtick_labels_product = [ligand.replace('HALO_S7P_', 'HALO_') for ligand in xtick_labels_product]
    plt.xticks(range(-1, len(dG_substrate_list)), xtick_labels_substrate, rotation=90)
    plt.sca(ax6)  # set axis current axis active
    plt.xticks(range(-1, len(dG_product_list)), xtick_labels_product, rotation=90)

    sns.boxplot(dG_substrate_list, ax=ax7, orient='v')
    sns.boxplot(dG_product_list, ax=ax8, orient='v')
    ax7.set_ylabel('dG (kcal/mol)')
    ax7.set_ylim(ymin, ymax)
    ax8.set_ylim(ymin, ymax)
    ax7.set_title(''.join([substrate, ' dG distribution']))
    ax8.set_title(''.join([product, ' dG distribution']))

    plt.tight_layout()
    plt.savefig(''.join([file_out, '.pdf']))
    plt.close()


def plot_dG_values_paper(data_df, file_out, substrate, product):
    """

    Plots all dG values in a scatter plot with error bars for the paper.

    :param data_df: dataframe with dG values for each MD simulation
    :param file_out: path+name for plot file
    :param substrate: substrate bound to the enzyme
    :param product: product bound to the enzyme
    :return: None
    """

    plot_timecourses_defs()

    dG_substrate_list = data_df[data_df['ligand'] == substrate]['mean'][:-1].values
    dG_std_substrate_list = data_df[data_df['ligand'] == substrate]['std'][:-1].values

    dG_product_list = data_df[data_df['ligand'] == product]['mean'][:-1].values
    dG_std_product_list = data_df[data_df['ligand'] == product]['std'][:-1].values

    dG_mean_substrate = data_df[data_df['ligand'] == substrate]['mean'].values[-1]
    dG_mean_product = data_df[data_df['ligand'] == product]['mean'].values[-1]

    xtick_labels_substrate = list(data_df[data_df['ligand'] == substrate].index.values[:-1])
    xtick_labels_substrate = [item[2:] for item in xtick_labels_substrate]
    xtick_labels_substrate.insert(0, '')
    #xtick_labels_substrate = ['$' + label.replace('_', '\_') +'$'  for label in xtick_labels_substrate]

    xtick_labels_product = list(data_df[data_df['ligand'] == product].index.values[:-1])
    xtick_labels_product = [item[2:] for item in xtick_labels_product]
    xtick_labels_product.insert(0, '')
    #xtick_labels_product = ['$' + label.replace('_', '\_') +'$' for label in xtick_labels_product]

    ymin = (np.min([np.min(dG_substrate_list-dG_std_substrate_list), np.min(dG_product_list-dG_std_product_list)])) - 5
    ymax = (np.max([np.max(dG_substrate_list+dG_std_substrate_list), np.max(dG_product_list+dG_std_product_list)])) + 5

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6, 6)) #, sharey=True)

    ax[0].errorbar(x=range(len(dG_substrate_list)), y=dG_substrate_list, yerr=dG_std_substrate_list, fmt='o', color='#377eb8', markersize=8)
    ax[0].axhline(y=dG_mean_substrate, color='black', linewidth=1.5)

    ax[1].errorbar(x=range(len(dG_product_list)), y=dG_product_list, yerr=dG_std_product_list, fmt='o',  color='#377eb8', markersize=8)
    ax[1].axhline(y=dG_mean_product, color='black', linewidth=1.5)

    ax[0].set_title(substrate)
    ax[1].set_title(product)

    ax[0].set_xlim(-1, len(dG_substrate_list))
    ax[1].set_xlim(-1, len(dG_product_list))
    ax[0].set_ylim(ymin, ymax)
    ax[1].set_ylim(ymin, ymax)
    #ax[0].set_xlabel('measurement')
    ax[1].set_xlabel('MD simulation')
    ax[0].set_ylabel('$\Delta G\,(kcal/mol)$')
    ax[1].set_ylabel('$\Delta G\,(kcal/mol)$')

    #print xtick_labels_substrate
    #print xtick_labels_product

    plt.sca(ax[0])  # set axis current axis active
    plt.xticks(range(-1, len(dG_substrate_list)), xtick_labels_substrate, rotation=90)
    plt.sca(ax[1])  # set axis current axis active
    plt.xticks(range(-1, len(dG_product_list)), xtick_labels_product, rotation=90)

    #GAPD and TALB
    #fig.subplots_adjust(top=0.95, bottom=0.23, left=0.16, right=0.99, hspace=1.1)
    # ENO_AB
    fig.subplots_adjust(top=0.95, bottom=0.17, left=0.16, right=0.99, hspace=0.7)
    #plt.tight_layout()
    plt.savefig(''.join([file_out, '.pdf']))
    plt.close()


def _get_mad(values_list):
    """

    Given a list of values, get the median absolute deviation: median(abs(value - median(values_list))

    :param values_list: list with values calculate the mad from
    :return mad: median absolute deviation of values_list
    """

    mad = np.median(abs(values_list - np.median(values_list)))

    return mad


def _analyze_and_plot_dG_values(data_df, file_out, substrate, product):

    dG_substrate_list = data_df[data_df['ligand'] == substrate]['mean'][:-1].values
    dG_std_substrate_list = data_df[data_df['ligand'] == substrate]['std'][:-1].values

    dG_product_list = data_df[data_df['ligand'] == product]['mean'][:-1].values
    dG_std_product_list = data_df[data_df['ligand'] == product]['std'][:-1].values

    dG_mean_substrate = data_df[data_df['ligand'] == substrate]['mean'][-1]
    dG_mean_product = data_df[data_df['ligand'] == product]['mean'][-1]

    xtick_labels_substrate = list(data_df[data_df['ligand'] == substrate].index.values[:-1])
    xtick_labels_substrate.insert(0, '')
    xtick_labels_product = list(data_df[data_df['ligand'] == product].index.values[:-1])
    xtick_labels_product.insert(0, '')

    ymin = (np.min([np.min(dG_substrate_list-dG_std_substrate_list), np.min(dG_product_list-dG_std_product_list)])) - 5
    ymax = (np.max([np.max(dG_substrate_list+dG_std_substrate_list), np.max(dG_product_list+dG_std_product_list)])) + 5

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10,7)) #, sharey=True)

    ax[0, 0].errorbar(x=range(len(dG_substrate_list)), y=dG_substrate_list, yerr=dG_std_substrate_list, fmt='o')
    ax[0, 0].axhline(y=dG_mean_substrate, color='black', linewidth=1.5)

    ax[0, 1].errorbar(x=range(len(dG_product_list)), y=dG_product_list, yerr=dG_std_product_list, fmt='o')
    ax[0, 1].axhline(y=dG_mean_product, color='black', linewidth=1.5)

    ax[0, 0].set_title(substrate)
    ax[0, 1].set_title(product)

    ax[0, 0].set_xlim(-1, len(dG_substrate_list))
    ax[0, 1].set_xlim(-1, len(dG_product_list))
    ax[0, 0].set_ylim(ymin, ymax)
    ax[0, 1].set_ylim(ymin, ymax)
    ax[0, 0].set_xlabel('measurement')
    ax[0, 1].set_xlabel('measurement')
    ax[0, 0].set_ylabel('dG (kcal/mol)')

    plt.sca(ax[0, 0])  # set axis current axis active
    plt.xticks(range(-1, len(dG_substrate_list)), xtick_labels_substrate, rotation=90)
    plt.sca(ax[0, 1])  # set axis current axis active
    plt.xticks(range(-1, len(dG_product_list)), xtick_labels_product, rotation=90)

    sns.boxplot(dG_substrate_list, ax=ax[1, 0], orient='v')
    sns.boxplot(dG_product_list, ax=ax[1, 1], orient='v')
    ax[1, 0].set_ylabel('dG (kcal/mol)')
    ax[1, 0].set_ylim(ymin, ymax)
    ax[1, 1].set_ylim(ymin, ymax)

    plt.tight_layout()
    plt.savefig(''.join([file_out, '.pdf']))
    plt.close()

    with open(''.join([file_out, '_tests.txt']), 'w') as f_out:
        f_out.write(str(spearmanr(range(len(dG_substrate_list)), dG_substrate_list-dG_mean_substrate)))
        f_out.write(str(spearmanr(range(len(dG_product_list)), dG_product_list-dG_mean_product)))
        f_out.write(str(shapiro(dG_substrate_list)))
        f_out.write(str(shapiro(dG_product_list)))


def get_final_results(form_list, ligand_list, cluster_list, root_dir, enzyme, start_frame, end_frame_orig, label,
                      substrate, product, file_out):
    """

    Gathers dG values for each MD simulation plus variance, std, median, median absolute deviation, etc.

    :param form_list: list with enzyme forms, e.g. 'APO', 'HALO'.
    :param ligand_list: list with ligands that were bound to enzyme. list is assumed to be the same for each enzyme form.
    :param cluster_list: list of lists with clusters for each enzyme form / ligand.
    :param root_dir: folder with all MMPBSA results
    :param enzyme: enzyme name
    :param start_frame: start frame to get dG values
    :param end_frame_orig: end frame to get dG values
    :param label: file label, e.g., 5ns
    :param substrate: substrate bound to enzyme
    :param product: product bound to enzyme
    :param file_out: path+name to file where dG will be stored
    :return: None
    """

    form_folder_list = ['_'.join([str(i+1), form]) for i, form in enumerate(form_list)]
    ligand_folder_list = ['_'.join([str(i+1), ligand]) for i, ligand in enumerate(ligand_list)]

    dG_dic = OrderedDict()
    for form_i in range(len(form_list)):
        for ligand_i in range(len(ligand_list)):

            for cluster_i, cluster in enumerate(cluster_list[form_i][ligand_i]):
                base_dir = ''.join([root_dir, form_folder_list[form_i], '/', ligand_folder_list[ligand_i], '/', cluster, '/'])
                file_name = ''.join([enzyme, '_WT_', form_list[form_i], '_', ligand_list[ligand_i], '_docked_', cluster])
                file_in = ''.join([base_dir, file_name, '.dat'])

                data_df = pd.read_csv(file_in, sep='\t', index_col=0)
                if enzyme == 'ENO_AB' and label == 'all' and cluster == 'cl000':
                    end_frame = 500
                else:
                    end_frame = end_frame_orig

                dG_values = data_df.ix[start_frame:end_frame, 0].values
                #print ligand_list[ligand_i], cluster, dG_values

                dG_dic['_'.join([form_list[form_i], ligand_list[ligand_i], cluster])] = [ligand_list[ligand_i], np.mean(dG_values), np.var(dG_values), np.std(dG_values),
                                                                                          np.var(dG_values), np.std(dG_values), np.median(dG_values), _get_mad(dG_values)]

    data_df = pd.DataFrame.from_dict(dG_dic, orient='index')
    data_df.columns = ['ligand', 'mean', 'var_propagated', 'std_propagated', 'var', 'std', 'median', 'mad']

    variance_of_average = get_variance_of_the_average(data_df[data_df['ligand'] == substrate]['var'])
    data_df.loc[len(data_df.index)] = [substrate,
                                       np.mean(data_df[data_df['ligand'] == substrate]['mean']),
                                       variance_of_average,
                                       np.sqrt(variance_of_average),
                                       np.var(data_df[data_df['ligand'] == substrate]['mean']),
                                       np.sqrt(np.var(data_df[data_df['ligand'] == substrate]['mean'])),
                                       np.median(data_df[data_df['ligand'] == substrate]['mean']),
                                       _get_mad(data_df[data_df['ligand'] == substrate]['mean'])]

    data_df.index.values[-1] = ''.join([substrate, ' sumup'])

    variance_of_average = get_variance_of_the_average(data_df[data_df['ligand'] == product]['var'])
    data_df.loc[len(data_df.index)] = [product,
                                       np.mean(data_df[data_df['ligand'] == product]['mean']),
                                       variance_of_average,
                                       np.sqrt(variance_of_average),
                                       np.var(data_df[data_df['ligand'] == product]['mean']),
                                       np.sqrt(np.var(data_df[data_df['ligand'] == product]['mean'])),
                                       np.median(data_df[data_df['ligand'] == product]['mean']),
                                       _get_mad(data_df[data_df['ligand'] == product]['mean'])]

    data_df.index.values[-1] = ''.join([product, ' sumup'])


    ddG = data_df.ix[-1, 1] - data_df.ix[-2, 1]
    ddG_variance_prop = data_df.ix[-1, 2] + data_df.ix[-2, 2]
    ddG_std_prop = np.sqrt(ddG_variance_prop)

    ddG_variance = data_df.ix[-1, 4] + data_df.ix[-2, 4]
    ddG_std = np.sqrt(ddG_variance)

    ddG_median = data_df.ix[-1, 6] - data_df.ix[-2, 6]
    ddG_mad = data_df.ix[-1, 7] + data_df.ix[-2, 7]

    data_df.loc[len(data_df.index)] = [''.join([product, '-', substrate]), ddG, ddG_variance_prop, ddG_std_prop, ddG_variance, ddG_std, ddG_median, ddG_mad]
    data_df.index.values[-1] = 'ddG'
    data_df.loc[len(data_df.index)] = [''.join([product, '-', substrate]), get_dKd_from_ddG(ddG), get_dKd_from_ddG(ddG-ddG_std_prop), get_dKd_from_ddG(ddG+ddG_std_prop),
                                       get_dKd_from_ddG(ddG-ddG_std), get_dKd_from_ddG(ddG+ddG_std), get_dKd_from_ddG(ddG_median), get_dKd_from_ddG(ddG_median+ddG_mad)]
    data_df.index.values[-1] = 'dKd'
    data_df['std_prop %'] = abs(data_df['std_propagated'] / data_df['mean'] * 100)
    data_df['std %'] = abs(data_df['std'] / data_df['mean'] * 100)

    file_out = '_'.join([file_out, label])
    data_df.to_csv(''.join([file_out, '.csv']), sep='\t')
    _analyze_and_plot_dG_values(data_df, file_out, substrate, product)


def check_prod_len_impact_on_ddG(file_in, prod_len_list, metric_list, metric_label_list):
    """

    Plot ddG/dKb/dG values for different MD simulation lengths.

    :param file_in: file with csv data
    :param prod_len_list: list with different MD simulations lengths in nanoseconds
    :param metric_list: list with metrics to calculate
    :param metric_label_list: list with metrics labels
    :return: None
    """

    metric_avg_dic = {metric: [] for metric in metric_list}
    metric_std_dic = {metric: [] for metric in metric_list}
    dKd_lb = []
    dKd_ub = []
    for prod_len in prod_len_list:
        data_df = pd.read_csv(''.join([file_in, '_', str(prod_len), 'ns.csv']), sep='\t', index_col=None, header=0)

        for i, metric in enumerate(metric_list):
            metric_avg_dic[metric].append(data_df[data_df.ix[:, 0] == metric]['mean'].values[0])
            metric_std_dic[metric].append(data_df[data_df.ix[:, 0] == metric]['std'].values[0])
            if metric == 'dKd':
                dKd_lb.append(data_df[data_df.ix[:, 0] == metric]['var'].values[0])
                dKd_ub.append(data_df[data_df.ix[:, 0] == metric]['std'].values[0])

    dKd_lb = np.array(metric_avg_dic['dKd']) - np.array(dKd_lb)
    dKd_ub = np.array(dKd_ub) - np.array(metric_avg_dic['dKd'])

    fig, ax = plt.subplots(nrows=len(metric_list), ncols=1, figsize=[6, 10])
    for i, metric in enumerate(metric_list):
        if metric == 'dKd':
            ax[i].errorbar(x=prod_len_list, y=metric_avg_dic[metric], yerr=[dKd_lb, dKd_ub], fmt='-o')
        else:
            ax[i].errorbar(x=prod_len_list, y=metric_avg_dic[metric], yerr=metric_std_dic[metric], fmt='-o')
        ax[i].set_xlabel('production length (ns)')
        ax[i].set_ylabel(''.join([metric_label_list[i], ' (kcal/mol)']))
        #if file_in.find('ENO') != -1 and metric == 'dKd':
        if metric == 'dKd':
            ax[i].set_yscale('log')

    plt.tight_layout()
    plt.savefig(''.join([file_in, '_evolution.pdf']))
    plt.close()


def gather_dG_plots(form_list, ligand_list, cluster_list, root_dir, enzyme, n_frames_per_form_per_ligand, file_out):
    """

    Gets dG values per frame for each enzyme form/ligand/cluster and plots them in a single document.

    :param form_list: list with enzyme forms considered, e.g., ['APO', 'HALO_S7P_remS7P']
    :param ligand_list: list with ligands considered per enzyme form  ['2PG', 'PEP']
    :param cluster_list: list of cluster per enzyme form per ligand
    :param root_dir: path to folder with all MMPBSA results
    :param enzyme: enzyme anme
    :param n_frames_per_form_per_ligand: how many frames used to calculate the dG per MD simulation
    :param file_out: path+name for plot file
    :return: None
    """

    form_folder_list = ['_'.join([str(i+1), form]) for i, form in enumerate(form_list)]
    ligand_folder_list = ['_'.join([str(i+1), ligand]) for i, ligand in enumerate(ligand_list)]

    n_rows = len(form_list)*n_frames_per_form_per_ligand
    fig, ax = plt.subplots(nrows=n_rows, ncols=len(ligand_list), figsize=[10, n_rows*3])
    min_dG = 10**10
    max_dG = -10**10
    for ligand_i in range(len(ligand_list)):
        for form_i in range(len(form_list)):
            for cluster_i, cluster in enumerate(cluster_list[form_i][ligand_i]):

                base_dir = ''.join([root_dir, form_folder_list[form_i], '/', ligand_folder_list[ligand_i], '/', cluster, '/'])
                file_name = ''.join([enzyme, '_WT_', form_list[form_i], '_', ligand_list[ligand_i], '_docked_', cluster])
                file_in = ''.join([base_dir, file_name, '.dat'])

                data_df = pd.read_csv(file_in, sep='\t', index_col=0)
                dG_values = data_df.ix[:, 0].values
                dG_values = dG_values[:-2]

                if max(dG_values) > max_dG:
                    max_dG = max(dG_values)
                if min(dG_values) < min_dG:
                    min_dG = min(dG_values)

                ax[form_i*n_frames_per_form_per_ligand + cluster_i][ligand_i].plot(dG_values)
                ax[form_i*n_frames_per_form_per_ligand + cluster_i][ligand_i].set_title(''.join([ligand_list[ligand_i], ', ', form_list[form_i], ', ', cluster]))

    for i in range(len(ax)):
        for j in range(len(ax[0])):
            ax[i][j].set_ylim(min_dG, max_dG)

    plt.tight_layout()
    plt.savefig(''.join([file_out, '.pdf']))
    plt.close()


def analyze_ENO(base_dir):

    enzyme = '1E9I'
    root_dir = ''.join([base_dir, '/ENO_1E9I/4_MMPBSA2016/'])
    form_list = ['APO']
    ligand_list = ['2PG', 'PEP']
    cluster_list = [[['cl000', 'cl000_2', 'cl000_3', 'cl000_4', 'cl100', 'cl100_2', 'cl100_3', 'cl100_4', 'cl200', 'cl200_2', 'cl200_3', 'cl200_4'],
                    ['cl000', 'cl000_2', 'cl000_3', 'cl000_4', 'cl100', 'cl100_2', 'cl100_3', 'cl100_4', 'cl110', 'cl110_2', 'cl110_3', 'cl110_4']]]
    file_out = ''.join([base_dir, '/ENO_1E9I/4_MMPBSA2016/ENO_ddG_sumup'])

    start_frame = 8
    end_frame_list = [35, 70, 105, 140, 175, 210]
    #end_frame_list = [140, 175, 210]
    label_list = ['5ns', '10ns', '15ns', '20ns', '25ns', '30ns']
    #label_list = ['20ns', '25ns', '30ns']

    for end_frame, label in zip(end_frame_list, label_list):
        get_final_results(form_list, ligand_list, cluster_list, root_dir, enzyme, start_frame, end_frame, label, ligand_list[0], ligand_list[1], file_out)

    data_df_rmsd = pd.read_csv(''.join([base_dir, '/ENO_1E9I/3_MD_post_dock2016/ENO_rmsd_sumup.csv']), sep='\t', index_col=0)
    data_df_dG = pd.read_csv(''.join([file_out, '_', label_list[-1], '.csv']), sep='\t', index_col=0)
    file_out_rmsd_dG = ''.join([base_dir, '/ENO_1E9I/4_MMPBSA2016/ENO_ddG_rmsd_sumup'])
    analyze_and_plot_dG_values_and_rmsd(data_df_rmsd, data_df_dG, file_out_rmsd_dG, ligand_list[0], ligand_list[1])

    file_in = file_out
    prod_len_list = [5, 10, 15, 20, 25, 30]
    #prod_len_list = [20, 25, 30]
    metric_list = ['ddG', 'dKd', '2PG sumup', 'PEP sumup']
    metric_label_list = ['ddG', 'dKd', 'dG 2PG', 'dG PEP']
    check_prod_len_impact_on_ddG(file_in, prod_len_list, metric_list, metric_label_list)

    file_out = ''.join([base_dir, '/ENO_1E9I/4_MMPBSA2016/ENO_all_dG_plots'])
    n_frames_per_form_per_ligand = 12
    gather_dG_plots(form_list, ligand_list, cluster_list, root_dir, enzyme, n_frames_per_form_per_ligand, file_out)


def analyze_eno_3h8a(base_dir):

    enzyme = 'ENO_AB'
    root_dir = ''.join([base_dir, '/ENO_AB/4_MMPBSA2016/'])
    form_list = ['MG']
    ligand_list = ['2PG', 'PEP']
    #cluster_list = [[['cl000', 'cl102', 'cl202', 'cl203'],
    #                ['cl000', 'cl302', 'cl312', 'cl402']]]

    cluster_list = [[['cl000', 'cl102', 'cl203'],
                    ['cl000', 'cl302', 'cl402']]]
    file_out = ''.join([base_dir, '/ENO_AB/4_MMPBSA2016/ENO_ddG_sumup'])

    start_frame = 8
    end_frame_list = [35, 70, 105, 139]
    label_list = ['5ns', '10ns', '15ns', '20ns']

    for end_frame, label in zip(end_frame_list, label_list):
        get_final_results(form_list, ligand_list, cluster_list, root_dir, enzyme, start_frame, end_frame, label, ligand_list[0], ligand_list[1], file_out)

    get_final_results(form_list, ligand_list, cluster_list, root_dir, enzyme, start_frame, 140, 'all', ligand_list[0], ligand_list[1], file_out)


    data_df_rmsd = pd.read_csv(''.join([base_dir, '/ENO_AB/3_MD_post_dock2016/ENO_rmsd_sumup.csv']), sep='\t', index_col=0)
    #data_df_dG = pd.read_csv(''.join([file_out, '_', label_list[-1], '.csv']), sep='\t', index_col=0)
    data_df_dG = pd.read_csv(''.join([file_out, '_all.csv']), sep='\t', index_col=0)
    file_out_rmsd_dG = ''.join([base_dir, '/ENO_AB/4_MMPBSA2016/ENO_ddG_rmsd_sumup'])
    analyze_and_plot_dG_values_and_rmsd(data_df_rmsd, data_df_dG, file_out_rmsd_dG, ligand_list[0], ligand_list[1])

    file_in = file_out
    prod_len_list = [5, 10, 15, 20]
    #prod_len_list = [20, 25, 30]
    metric_list = ['ddG', 'dKd', '2PG sumup', 'PEP sumup']
    metric_label_list = ['ddG', 'dKd', 'dG 2PG', 'dG PEP']
    check_prod_len_impact_on_ddG(file_in, prod_len_list, metric_list, metric_label_list)

    file_out = ''.join([base_dir, '/ENO_AB/4_MMPBSA2016/ENO_all_dG_plots'])
    n_frames_per_form_per_ligand = 4
    gather_dG_plots(form_list, ligand_list, cluster_list, root_dir, enzyme, n_frames_per_form_per_ligand, file_out)


def analyze_GAPD(base_dir):

    enzyme = 'GAPDH'
    root_dir = ''.join([base_dir, '/GAPDH/4_MMPBSA2016/'])
    file_out = ''.join([base_dir, '/GAPDH/4_MMPBSA2016/GAPD_ddG_sumup'])
    form_list = ['NAD', 'NAD_remG3P']
    ligand_list = ['DPG']
    cluster_list = [[['cl100', 'cl100_2', 'cl100_3', 'cl100_4', 'cl100_5',
                        'cl101', 'cl101_2', 'cl101_3', 'cl101_4', 'cl101_5',
                        'cl120', 'cl120_2', 'cl120_3', 'cl120_4', 'cl120_5',
                        'cl220', 'cl220_2', 'cl220_3', 'cl220_4', 'cl220_5']]]

    start_frame = 8
    end_frame_list = [70]
    label_list = ['1ns']

    for end_frame, label in zip(end_frame_list, label_list):
        get_final_results(form_list, ligand_list, cluster_list, root_dir, enzyme, start_frame, end_frame, label, ligand_list[0], ligand_list[1], file_out)

    #data_df_rmsd = pd.read_csv(''.join([base_dir, '/GAPDH/3_MD_post_dock2016/GAPD_rmsd_sumup.csv']), sep='\t', index_col=0)
    #data_df_dG = pd.read_csv(''.join([file_out, '_', label_list[-1], '.csv']), sep='\t', index_col=0)
    #file_out_rmsd_dG = ''.join([base_dir, '/GAPDH/4_MMPBSA2016/GAPD_ddG_rmsd_sumup'])
    #analyze_and_plot_dG_values_and_rmsd(data_df_rmsd, data_df_dG, file_out_rmsd_dG, ligand_list[0], ligand_list[1])


    #file_out = ''.join([base_dir, '/GAPDH/4_MMPBSA2016/GAPD_all_dG_plots'])
    #n_frames_per_form_per_ligand = 6
    #gather_dG_plots(form_list, ligand_list, cluster_list, root_dir, enzyme, n_frames_per_form_per_ligand, file_out)


def analyze_TALB(base_dir):

    enzyme = 'TALB'
    root_dir = ''.join([base_dir, '/TalB/4_MMPBSA2016/'])
    file_out = ''.join([base_dir, '/TalB/4_MMPBSA2016/TALB_ddG_sumup'])
    form_list = ['APO', 'HALO_S7P_remS7P']
    ligand_list = ['S7P_linear', 'F6P']
    cluster_list = [[['cl101', 'cl101_2', 'cl102', 'cl102_2', 'cl212', 'cl212_2'],
                     ['cl202', 'cl202_2', 'cl212', 'cl212_2', 'cl321', 'cl321_2']],
                    [['cl101', 'cl101_2', 'cl302', 'cl302_2', 'cl322', 'cl322_2'],
                     ['cl101', 'cl101_2', 'cl201', 'cl201_2', 'cl311', 'cl311_2']]]

    start_frame = 8
    end_frame_list = [35, 70, 105, 139]
    label_list = ['5ns', '10ns', '15ns', '20ns']

    for end_frame, label in zip(end_frame_list, label_list):
        get_final_results(form_list, ligand_list, cluster_list, root_dir, enzyme, start_frame, end_frame, label, ligand_list[0], ligand_list[1], file_out)

    data_df_rmsd = pd.read_csv(''.join([base_dir, '/TalB/3_MD_post_dock2016/TALB_rmsd_sumup.csv']), sep='\t', index_col=0)
    data_df_dG = pd.read_csv(''.join([file_out, '_', label_list[-1], '.csv']), sep='\t', index_col=0)
    file_out_rmsd_dG = ''.join([base_dir, '/TalB/4_MMPBSA2016/TALB_ddG_rmsd_sumup'])
    analyze_and_plot_dG_values_and_rmsd(data_df_rmsd, data_df_dG, file_out_rmsd_dG, ligand_list[0], ligand_list[1])

    file_in = file_out
    prod_len_list = [5, 10, 15, 20]
    metric_list = ['ddG', 'dKd', 'S7P_linear sumup', 'F6P sumup']
    metric_label_list = ['ddG', 'dKd', 'dG S7P', 'dG F6P']
    check_prod_len_impact_on_ddG(file_in, prod_len_list, metric_list, metric_label_list)

    file_out = ''.join([base_dir, '/TalB/4_MMPBSA2016/TALB_all_dG_plots'])
    n_frames_per_form_per_ligand = 6
    gather_dG_plots(form_list, ligand_list, cluster_list, root_dir, enzyme, n_frames_per_form_per_ligand, file_out)


def main():
    base_dir = '/home/mrama/Dropbox/PhD_stuff/Projects/MD/eMASS-MD/MD_data'

    analyze_ENO(base_dir)

    analyze_eno_3h8a(base_dir)

    analyze_GAPD(base_dir)

    analyze_TALB(base_dir)


if __name__ == '__main__':
    main()