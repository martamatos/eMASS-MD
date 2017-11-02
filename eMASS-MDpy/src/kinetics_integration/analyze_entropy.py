import math
import re
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from src.kinetics_integration.binning import bin_keqs
from src.kinetics_integration.plots_definitions import entropy_plot_single_defs
from src.kinetics_integration.utils import import_rateconstants


def calculate_set_entropy(data_df, n_samples):
    n_patterns = len(data_df.index)
    entropy = 0
    for i in range(n_patterns):
        p_xi = data_df.ix[i, 0] / float(n_samples)
        # print p_xi
        entropy += p_xi * math.log10(p_xi)
    entropy = -entropy
    return entropy


def _plot_entropy(entropy_list, entropy_label, fixed_size, bin_width_list, file_out_base, column_labels, y_lims, limit):
    print(entropy_list)

    n_rows = len(bin_width_list) / 2 if len(bin_width_list) % 2 == 0 else len(bin_width_list) / 2 + 1
    n_cols = 2
    fig_size = (10, 4 * n_rows)
    fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=fig_size)

    for i in range(n_rows):
        for j in range(n_cols):
            ind = i * n_cols + j
            if ind >= len(entropy_list):
                break
            #print entropy_list[ind].columns
            if column_labels:
                entropy_list[ind].columns = column_labels

            if n_rows == 1:
                ax_ij = ax[j]
            else:
                ax_ij = ax[i, j]
            sns.boxplot(entropy_list[ind], ax=ax_ij)
            plt.sca(ax_ij)  # set axis current axis active
            plt.xticks(rotation=90)

            if re.search('rel', entropy_label):
                ax_ij.set_ylim(y_lims)
            else:
                y_max = entropy_list[ind].max().max()
                ax_ij.set_ylim(0, y_max)
            ax_ij.set_xlabel('data removed from fitting')
            # ax[i, j].set_xlabel('data used for fitting')
            ax_ij.set_ylabel('entropy')

            title = ''.join(['bin width (orders of magnitude): ', str(bin_width_list[ind])])
            ax_ij.set_title(title)

    plt.tight_layout()
    if fixed_size:
        file_out = ''.join([file_out_base, str(limit), '_', entropy_label, '_entropy_fixed_size.'])
    else:
        file_out = ''.join([file_out_base, str(limit), '_', entropy_label, '_entropy_fixed_width.'])
    plt.savefig(''.join([file_out, 'pdf']))
    plt.savefig(''.join([file_out, 'png']))
    plt.close()


def _plot_entropy_single(entropy_list, entropy_label, fixed_size, bin_width_list, file_out_base, column_labels, y_lims,
                         x_label, color_list):

    entropy_plot_single_defs()

    fig, ax = plt.subplots(figsize=(6.5, 4.2))

    entropy = entropy_list[0]  # works because entropy_list only has one entry, otherwise _plot_entropy() would be used

    entropy.columns = column_labels

    sns.boxplot(data=entropy, ax=ax, palette=color_list, orient='v')
    plt.xticks(rotation=0)
    ax.set_ylim(y_lims)

    ref_median = np.median(entropy['None'].dropna())
    dKd_median = np.median(entropy['$\Delta K_d$'].dropna())

    plt.plot((-1, 9), (ref_median, ref_median), linestyle='--', color='#e41a1c')
    plt.plot((-1, 9), (dKd_median, dKd_median), linestyle='--', color='#377eb8')

    ax.set_xlabel(x_label)
    ax.set_ylabel('Relative Shannon entropy')

    plt.tight_layout()
    if fixed_size:
        file_out = ''.join([file_out_base, entropy_label, '_entropy_fixed_size.'])
    else:
        file_out = ''.join(
            [file_out_base, entropy_label, '_entropy_fixed_width_', str(bin_width_list[0]),
             '.'])
    plt.savefig(''.join([file_out, 'pdf']))
    plt.savefig(''.join([file_out, 'png']), dpi=100)
    plt.close()


def calculate_enzyme_entropy(file_in_base, file_out_base, n_model_ensembles, model_types, column_order, n_variables,
                             ssd_threshold, Keq, convert_to_ratios, limit, bin_width_list, column_labels, y_lims,
                             x_label, color_list=None, plot_entropy_only=False, fixed_size=False,
                             n_samples_per_bin=None):
    """

    Given a set of model ensembles, calculates and plots the relative Shannon entropy for each model ensemble.

    :param file_in_base: path and base name for the input file
    :param file_out_base: path and base name for the output file
    :param n_model_ensembles: integer specifying how many model ensembles are to be analyzed
    :param model_types: a list with all model types, e.g. ['all', 'dKd']
    :param column_order: list with the order the columns in the file to be imported should follow
    :param n_variables: how many rate or elementary equilibrium constants in each model
    :param ssd_threshold: threshold for sum of squared deviations, models where ssd > ssd_threshold are discarded
    :param Keq: equilibrium constant value
    :param convert_to_ratios: boolean for whether the rate constants should be converted into elementary equilibrium
                              constants or not
    :param limit: integer specifying how many models should be considered
    :param bin_width_list: list with different bin widths (in order of magnitude)
    :param column_labels: list with labels for each column in the file to be imported
    :param y_lims: tuple or list with limits for the y-axis
    :param x_label: x-axis label
    :param color_list: list with colors to use on the plot
    :param plot_entropy_only: boolean specifying
    :param fixed_size: boolean specifying whether the bins should have a fixed size, i.e. a fixed number of samples per
                       bin.
    :param n_samples_per_bin: how many samples per bin if fixed_size=True
    :return:
    """

    n_models = OrderedDict()
    rel_entropy_list = []
    rel_model_entropy = OrderedDict()

    if not plot_entropy_only:

        for i, bin_width in enumerate(bin_width_list):

            #print 'bin width', bin_width
            if file_out_base.find('ratios') != -1:
                n_bins = 15 / bin_width if 15 % bin_width == 0 else 15 / bin_width + 1
                n_bins *= 2
            else:
                n_bins = 15 / bin_width if 15 % bin_width == 0 else 15 / bin_width + 1

            for model_type in model_types:
                rel_model_entropy[model_type] = []
                n_models[model_type] = []

                #print 'model_type', model_type
                for ensemble_i in range(n_model_ensembles):

                    file_in = ''.join([file_in_base, model_type, '_', str(ensemble_i), '.csv'])

                    if re.search('Keq', model_type):
                        Keq_local = None
                    else:
                        Keq_local = Keq
                    try:
                        data_df, fitness = import_rateconstants(file_in, ssd_threshold, Keq_local, column_order,
                                                                convert_to_ratios=convert_to_ratios, limit=limit)
                    except IOError:
                        print file_in + ' does not exist'
                        n_models[model_type].append(0)
                        rel_model_entropy[model_type].append(None)

                        continue

                    n_samples = len(data_df.index)

                    n_models[model_type].append(n_samples)
                    if fixed_size:
                        file_out = ''.join(
                            [file_out_base, 'fixed_size_', str(n_samples_per_bin), '_', model_type, '_',
                             str(ensemble_i)])
                        pattern_count_df, pattern_df_by_model = bin_keqs(data_df, fixed_size=fixed_size,
                                                                         n_samples_per_bin=None, file_out=file_out)
                    elif bin_width:
                        pattern_count_df, pattern_df_by_model = bin_keqs(data_df, bin_width=bin_width)
                    else:
                        pattern_count_df, pattern_df_by_model = bin_keqs(data_df)

                    #print pattern_count_df
                    #exit()
                    assert n_samples == pattern_count_df.ix[:, 0].sum()
                    entropy = calculate_set_entropy(pattern_count_df, n_samples)
                    max_entropy = math.log10(n_bins ** n_variables)

                    rel_model_entropy[model_type].append(entropy / max_entropy)

                    #if limit < 100:
                    #    assert n_samples == limit

            n_models_df = pd.DataFrame().from_dict(n_models)
            n_models_df = n_models_df.transpose()
            n_models_df.to_csv(''.join([file_out_base, 'n_models.csv']), sep='\t')

            rel_entropy_df = pd.DataFrame().from_dict(rel_model_entropy)

            rel_entropy_df.to_csv(
                ''.join([file_out_base, str(limit), '_rel_entropy_fixed_width_', str(bin_width), '.csv']), sep='\t')

            rel_entropy_list.append(rel_entropy_df)

    else:
        for i, bin_width in enumerate(bin_width_list):
            rel_entropy_df = pd.read_csv(
                ''.join([file_out_base, str(limit), '_rel_entropy_fixed_width_', str(bin_width), '.csv']), sep='\t',
                index_col=0)

            rel_entropy_list.append(rel_entropy_df)
    #print(rel_entropy_list)

    if len(bin_width_list) > 1:
        entropy_label = ''.join(['rel_', str(limit), '_models'])
        _plot_entropy(rel_entropy_list, entropy_label, fixed_size, bin_width_list, file_out_base, column_labels, y_lims,
                      limit)
    if len(bin_width_list) == 1:
        entropy_label = ''.join(['rel_', str(limit), '_models'])
        _plot_entropy_single(rel_entropy_list, entropy_label, fixed_size, bin_width_list, file_out_base, column_labels,
                             y_lims, x_label, color_list)
