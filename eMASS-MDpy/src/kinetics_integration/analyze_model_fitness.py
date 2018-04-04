import re
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from src.kinetics_integration.plots_definitions import entropy_plot_single_defs
from src.kinetics_integration.utils import import_rateconstants


def plot_n_good_models_dist(file_in_base, file_out, model_type_list, n_ensembles, ssd_threshold, Keq,
                            column_order, column_labels, x_label, scale_data, convert_to_ratios, color_list, filter):

    """
    Given a file_in_base and model_type, it imports all the model ensembles for that model type, takes the fitness
    and plots it for all models in the n model ensembles of each model type.

    :param file_in_base: path and base name for the input file
    :param file_out: path and base name for the output file
    :param model_type_list: a list with all model types, e.g. ['all', 'dKd']
    :param n_ensembles: how many model ensembles to be analyzed
    :param ssd_threshold: threshold for sum of squared deviations
    :param Keq: equilibrium constant value
    :param column_order: list with the order the columns in the file to be imported should follow
    :param column_labels: list with labels for each column in the file to be imported
    :param x_label: x-axis label
    :param scale_data: boolean for whether or not the data should be scaled - not used anymore
    :param convert_to_ratios: boolean for whether the rate constants should be converted into elementary equilibrium
                              constants or not
    :param color_list: list with colors to use on the plot
    :param filter: boolean defining whether or not the models should be filtered by ssd
    :return: None
    """

    n_models_dic = OrderedDict((key, []) for key in model_type_list)

    n_model_types = len(model_type_list)

    entropy_plot_single_defs()

    n_rows = n_model_types // 3 if n_model_types % 3 == 0 else n_model_types // 3 + 1
    fig, ax = plt.subplots(figsize=[6, 4])

    for i, model_type in enumerate(model_type_list):

        n_models_list = []

        if re.search('Keq', model_type):
            Keq_local = None
        else:
            Keq_local = Keq

        for ensemble_i in range(n_ensembles):

            file_in = ''.join([file_in_base, model_type, '_', str(ensemble_i), '.csv'])
            try:
                data_df, fitness = import_rateconstants(file_in, ssd_threshold, Keq_local, column_order,
                                                        convert_to_ratios=convert_to_ratios, filter=filter)
            except IOError:
                print file_in + ' does not exist'
                continue

            n_models_list.append(len(data_df.index))

        n_models_dic[model_type] = n_models_list

    sns.boxplot(data=[n_models_dic[key] for key in model_type_list], ax=ax, palette=color_list)
    ax.set_ylim(0, 110)
    ax.set_xticklabels(column_labels, rotation=0)
    ax.set_xlabel(x_label)
    ax.set_ylabel('number of models with ssld$<0.1$')
    ax.grid(True)

    plt.tight_layout()
    plt.savefig(''.join([file_out, '.pdf']), dpi=300)
    plt.savefig(''.join([file_out, '.png']), dpi=100)
    plt.close()


def plot_fitness_dist(file_in_base, file_out, model_type_list, n_ensembles, ssd_threshold, Keq, column_order,
                      column_labels, x_label, convert_to_ratios, color_list, filter):

    """
    Given a file_in_base and model_type, it imports all the model ensembles for that model type, takes the fitness
    and plots it for all models in the n model ensembles of each model type.

    :param file_in_base: path and base name for the input file
    :param file_out: path and base name for the output file
    :param model_type_list: a list with all model types, e.g. ['all', 'dKd']
    :param n_ensembles: how many model ensembles to be analyzed
    :param ssd_threshold: threshold for sum of squared deviations
    :param Keq: equilibrium constant value
    :param column_order: list with the order the columns in the file to be imported should follow
    :param column_labels: list with labels for each column in the file to be imported
    :param x_label: x-axis label
    :param convert_to_ratios: boolean for whether the rate constants should be converted into elementary equilibrium
                              constants or not
    :param color_list: list with colors to use on the plot
    :param filter: boolean defining whether or not the models should be filtered by ssd
    :return: None
    """

    fitness_dic = OrderedDict((key, []) for key in model_type_list)

    n_model_types = len(model_type_list)

    entropy_plot_single_defs()

    n_rows = n_model_types // 3 if n_model_types % 3 == 0 else n_model_types // 3 + 1
    fig, ax = plt.subplots(figsize=[6, 4])

    for i, model_type in enumerate(model_type_list):

        fitness_list = []
        for ensemble_i in range(n_ensembles):

            file_in = ''.join([file_in_base, model_type, '_', str(ensemble_i), '.csv'])
            try:
                data_df, fitness = import_rateconstants(file_in, ssd_threshold, Keq, column_order,
                                                        convert_to_ratios=convert_to_ratios, filter=filter)
            except IOError:
                print file_in + ' does not exist'
                continue
            fitness_list.append(fitness)

        fitness_all = [item for sublist in fitness_list for item in sublist]
        fitness_dic[model_type] = fitness_all

    sns.boxplot(data=[np.log10(fitness_dic[key]) for key in model_type_list], ax=ax, palette=color_list)
    #sns.boxplot(data=[fitness_dic[key] for key in model_type_list], ax=ax)
    ax.set_xticklabels(column_labels, rotation=0)
    ax.set_xlabel(x_label)
    ax.set_ylabel('ssld order of magnitude ($10^y)$')
    ax.grid(True)
    plt.plot((-1, 9), (np.log10(0.1), np.log10(0.1)), color='black')

    plt.tight_layout()
    #print file_out
    plt.savefig(''.join([file_out, '.pdf']), dpi=300)
    plt.savefig(''.join([file_out, '.png']), dpi=100)
    plt.close()

