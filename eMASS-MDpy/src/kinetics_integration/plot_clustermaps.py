import matplotlib.pyplot as plt
import seaborn as sns

from src.kinetics_integration.plots_definitions import clustermap_defs
from src.kinetics_integration.utils import import_rateconstants


def _cluster_map_general(data_df, file_out, v_min, v_max, fig_size=(4, 6)):

    clustermap_defs()

    fig, ax = plt.subplots(1, 1)
    res = sns.clustermap(data=data_df, metric='euclidean', col_cluster=False, vmin=v_min, vmax=v_max, figsize=fig_size,
                         yticklabels=False, cmap='Reds', cbar_kws={'ticks': [-6, -3, 0, 3, 6, 9]})

    plt.tick_params(labelsize=18)
    plt.savefig(''.join([file_out, '_clustermap.pdf']))
    plt.savefig(''.join([file_out, '_clustermap.png']))

    plt.close()


def cluster_map_param_inf(file_in_base, file_out_base, ssd_threshold, convert_to_ratios, column_order, vmin, vmax,
                          model_types, n_ensembles, fig_size=(4, 6), Keq=None):

    """

    Plots a clustermap for each model ensemble for each model type, where different rate constant sets (or elementary
    equilibrium constant sets) are clustered together.

    :param file_in_base: path and base name for the input file
    :param file_out: path and base name for the output file
    :param ssd_threshold: threshold for sum of squared deviations
    :param convert_to_ratios: boolean for whether the rate constants should be converted into elementary equilibrium
                              constants or not
    :param column_order: list with the order the columns in the file to be imported should follow
    :param vmin: minimum value to anchor the colormap in the clustermap
    :param vmax: maximum value to anchor the colormap in the clustermap
    :param model_types: a list with all model types, e.g. ['all', 'dKd']
    :param n_ensembles: how many model ensembles to be analyzed
    :param Keq: equilibrium constant value
    :return: None

    """

    for ensemble_i in range(n_ensembles):

        for model_type in model_types:

            # print model_type
            file_in = ''.join([file_in_base, model_type, '_', str(ensemble_i), '.csv'])

            data_df, fitness = import_rateconstants(file_in, ssd_threshold, Keq, column_order,
                                                    convert_to_ratios=convert_to_ratios)

            file_out = ''.join([file_out_base, model_type, '_', str(ensemble_i)])
            _cluster_map_general(data_df, file_out, vmin, vmax, fig_size)

