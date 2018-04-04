from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from src.kinetics_integration.plots_definitions import get_elementary_keq_range_scatter_defs
from src.kinetics_integration.utils import import_rateconstants


def _get_keq_range(data_df, ekeq_range_dic, model_type):

    for i in range(len(data_df.columns)):
        q25 = np.percentile(data_df.ix[:, i], 25)
        q75 = np.percentile(data_df.ix[:, i], 75)

        ekeq_range_dic[model_type][i].append(q75-q25)


def _get_mad(ekeq_range_dic, key, y_data, n_keqs):
    mad_list = []
    for keq_i in range(n_keqs):
        mad_list.append(np.median(np.abs(y_data[keq_i] - ekeq_range_dic[key][keq_i])))

    return mad_list


def _plot_keq_iqr(ekeq_range_dic, n_model_types, n_keqs, column_labels, file_out, x_labels, color_list,
                  plot_type='boxplot', zorder_list=None):

    get_elementary_keq_range_scatter_defs()

    if plot_type == 'scatter':
        if not zorder_list:
            zorder_list = range(n_model_types)
        #print zorder_list
        fig, ax = plt.subplots(figsize=[6, 4.2])

        for i, key in enumerate(ekeq_range_dic.keys()):
            #print i,  key

            y_data = [np.median(ekeq_range_dic[key][keq_i]) for keq_i in range(n_keqs)]
            mad_list = _get_mad(ekeq_range_dic, key, y_data, n_keqs)
            #print y_data
            #print mad_list
            y_error_lb = [y_data[keq_i] - mad_list[keq_i] for keq_i in range(n_keqs)]
            y_error_ub = [mad_list[keq_i] + y_data[keq_i] for keq_i in range(n_keqs)]
            #print y_error_lb
            #print y_error_ub

            (_, caps, eline) = plt.errorbar(x=range(n_keqs), y=y_data, yerr=mad_list, fmt='-o',
                                            label=column_labels[i], markersize=8, linewidth=1, elinewidth=0.2,
                                            color=color_list[i], capthick=1, zorder=zorder_list[i])

            eline[0].set_linestyle('--')

            ax.set_xlim(-0.5, n_keqs-0.5)
            ax.set_ylim(0, 12)
        #print len(x_labels)
        plt.xticks(np.arange(0, len(x_labels), 1))
        #print x_labels
        ax.set_xticklabels(x_labels, rotation=0)

        ax.set_xlabel('Rate constants')
        ax.set_ylabel('Order of magnitude of $k$ range')

        ax.legend([plt.Line2D([], [], marker='o', color=color_list[i], markersize=9) for i in range(len(column_labels))],
                   column_labels,
                   numpoints=1,
                   ncol=4,
                   loc='upper center',
                   bbox_to_anchor=(0.5, 1.25))

        fig.subplots_adjust(top=0.85, bottom=0.2, left=0.12, right=0.99, hspace=0.7)
        plt.savefig(''.join([file_out, '_scatter.pdf']))
        plt.savefig(''.join([file_out, '_scatter.png']))
        plt.close()

    if plot_type == 'boxplot':

        n_rows = n_model_types // 3 if n_model_types % 3 == 0 else n_model_types // 3 + 1
        fig, ax = plt.subplots(nrows=n_rows, ncols=3, figsize=[12, n_rows * 3])

        if n_rows == 1:
            for i, key in enumerate(ekeq_range_dic.keys()):
                sns.boxplot(data=ekeq_range_dic[key], ax=ax[i])
                ax[i].set_title(column_labels[i])
                ax[i].set_ylim(0, 30)
        else:
            for i, key in enumerate(ekeq_range_dic.keys()):
                ind_j = i % 3
                ind_i = i // 3
                #print i, ind_i, ind_j
                sns.boxplot(data=ekeq_range_dic[key], ax=ax[ind_i, ind_j])
                ax[ind_i, ind_j].set_title(column_labels[i])
                ax[ind_i, ind_j].set_ylim(0, 30)

        plt.tight_layout()
        plt.savefig(''.join([file_out, '_boxplot.pdf']))
        plt.savefig(''.join([file_out, '_boxplot.png']))
        plt.close()


def calculate_elementary_keqs_iqr(file_in_base, file_out_base, model_type_list, n_model_ensembles, ssd_threshold, Keq,
                                  column_order, convert_to_ratios, column_labels, x_labels, color_list,
                                  zorder_list=None):
    """
    Calculates the inter quartile range for each rate constant or elementary equilibrium constant
    (if convert_to_ratios=True) and plots as both a boxplot or a scatter plot.


    :param file_in_base: path and base name for the input file
    :param file_out_base: path and base name for the output file
    :param model_type_list: a list with all model types, e.g. ['all', 'dKd']
    :param n_model_ensembles: integer specifying how many model ensembles are to be analyzed
    :param ssd_threshold: threshold for sum of squared deviations
    :param Keq: equilibrium constant value
    :param column_order: list with the order the columns in the file to be imported should follow
    :param convert_to_ratios: boolean for whether the rate constants should be converted into elementary equilibrium
                              constants or not
    :param column_labels: list with labels for each column in the file to be imported
    :param x_labels: x-axis label
    :param color_list: list with colors to use on the plot
    :param zorder_list: order in which lines are plotted in the graph
    :return: None
    """

    ekeq_range_dic = OrderedDict((model_type, [[] for i in range(len(column_order))]) for model_type in model_type_list)

    for model_type in model_type_list:
        for ensemble_i in range( n_model_ensembles):
            file_in = ''.join([file_in_base, model_type, '_', str(ensemble_i), '.csv'])
            #print model_type, model_i
            try:
                data_df, fitness = import_rateconstants(file_in, ssd_threshold, Keq, column_order,
                                                        convert_to_ratios=convert_to_ratios)
            except IOError:
                print file_in + ' does not exist'
                continue

            _get_keq_range(data_df, ekeq_range_dic, model_type)

    _plot_keq_iqr(ekeq_range_dic, len(model_type_list), len(data_df.columns), column_labels, file_out_base, x_labels,
                  color_list, plot_type='boxplot', zorder_list=zorder_list)
    _plot_keq_iqr(ekeq_range_dic, len(model_type_list), len(data_df.columns), column_labels, file_out_base, x_labels,
                  color_list, plot_type='scatter', zorder_list=zorder_list)



