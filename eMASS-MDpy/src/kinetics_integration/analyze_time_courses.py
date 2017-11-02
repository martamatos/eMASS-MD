import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter

from src.kinetics_integration.plots_definitions import plot_timecourses_defs, entropy_plot_single_defs
from src.kinetics_integration.utils import import_data_timecourses

COLOR_LIST = ['#377eb8', '#e41a1c',  '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']


def plot_timecourses(data_df_nodKd, data_df_dKd, file_out_base, species_list, ensemble_i, y_lims, Keq='',
                     mets_order=None, legend_label=None):
    """
    Given dataframes with and without dKd, plots data in a figure with 2 subplots.

    :param data_df_nodKd: dataframe with data without dKd
    :param data_df_dKd: dataframe with data with dKd
    :param file_out_base: path and base name for the output file
    :param species_list: list with species to be plotted: metabolites, enzymes, or reactions fluxes
    :param ensemble_i: which model ensemble is being plotted
    :param y_lims: list or tuple with limits for y-axis
    :param Keq: equilibrium constant value
    :param mets_order: list with order for metabolite/enzyme concentrations and reactions fluxes
    :param legend_label: label for legend
    :return: None
    """

    plot_timecourses_defs()

    if not mets_order:
        mets_order = species_list

    if not legend_label and mets_order:
        legend_label = mets_order
    elif not legend_label and not mets_order:
        legend_label = species_list

    n_species = len(species_list)
    fig, ax_list = plt.subplots(n_species, 2,  figsize=(6.5, 5), sharex=True, sharey=False)

    for i, met in enumerate(mets_order):
        ax_list[i, 0].plot(data_df_nodKd['t'].values, data_df_nodKd[met].values, c=COLOR_LIST[i])
        if y_lims:
            ax_list[i, 0].set_ylim(y_lims[i])

        ax_list[i, 1].plot(data_df_dKd['t'].values, data_df_dKd[met].values, c=COLOR_LIST[i])
        if y_lims:
            ax_list[i, 1].set_ylim(y_lims[i])
        ax_list[i, 1].set_yticklabels('')

    plt.sca(ax_list[0, 0])
    fig.gca().get_yaxis().set_major_formatter(ScalarFormatter(useOffset=False))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.sca(ax_list[1, 0])
    fig.gca().get_yaxis().set_major_formatter(ScalarFormatter(useOffset=False))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.xlabel("Time (s)")
    plt.ylabel("Concentration (mol/L)")

    fig.legend(
        [plt.Line2D((0, 1), (0, 0), color=COLOR_LIST[i], linestyle='-', linewidth=2) for i in range(len(species_list))],
        legend_label,
        numpoints=1,
        ncol=2,
        loc='center',
        bbox_to_anchor=(0.55, 0.95))

    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.12, right=0.96, wspace=0.16, hspace=0.2)
    plt.savefig(''.join([file_out_base, '_', str(ensemble_i), '_', str(Keq), '.pdf']))
    plt.close()


def _subplot_filtered(ax, data_df, species_list, x_lims, y_lims, title):
    for i, met in enumerate(species_list):
        ax.plot(data_df['t'].values, data_df[met].values, c=COLOR_LIST[i])

    ax.set_xlim(x_lims)
    ax.set_ylim(y_lims)
    ax.set_xlabel('Time (s)')
    # ax.set_ylabel('Flux (mol/L.s)')
    ax.set_ylabel('Concentration (mol/L)')
    ax.set_title(title)


def get_species_median(data_df, dic_dKd, dic_nodKd, offset=1):

    for species in data_df.columns.values[:-offset]:
        dic_dKd[species].append(np.median(data_df[data_df['model_type'] == 'dKd'][species]))
        dic_nodKd[species].append(np.median(data_df[data_df['model_type'] == 'nodKd'][species]))


def _plot_single(data_df, model_types, species_list, colorlist_swapped, file_out_base, y_axis_label, yticks_label, plot_spacings, fig_size=(6,5)):

    entropy_plot_single_defs()

    fig, ax = plt.subplots(figsize=fig_size)

    flierprops = dict(linestyle='none')

    sns.boxplot(x='species', y='concentration', hue='model_type', data=data_df,
                ax=ax, palette=colorlist_swapped, flierprops=flierprops)
    ax.legend_.remove()

    fig.gca().get_yaxis().set_major_formatter(ScalarFormatter(useOffset=False))
    fig.gca().get_yaxis().set_major_formatter(mtick.FormatStrFormatter('%.2f'))

    if yticks_label:
        ax.set_xticklabels(yticks_label, rotation=60)

    ax.set_xlabel('')
    ax.set_ylabel(y_axis_label)

    ax.yaxis.grid(True)

    fig.legend([plt.Line2D((0, 1), (0, 0), color=colorlist_swapped[i], linestyle='-', linewidth=2) for i in range(len(model_types))],
               ['no $\Delta K_b $', '$\Delta K_b$'],
               numpoints=1,
               ncol=2,
               loc='center',
               bbox_to_anchor=(0.55, 0.95))

    fig.subplots_adjust(top=plot_spacings[0], bottom=plot_spacings[1], left=plot_spacings[2], right=plot_spacings[3], wspace=0.4, hspace=0.5)
    #plt.tight_layout()
    #plt.xticks(rotation=70)
    plt.savefig(''.join([file_out_base, '_single.pdf']))
    plt.close()

    return fig, ax


def _plot_subplots(data_df, model_types, species_list, colorlist_swapped, file_out_base, y_label):

    n_species = len(species_list)
    n_rows = n_species // 3 if n_species % 3 == 0 else (n_species // 3)+1
    fig, ax = plt.subplots(n_rows, 3, figsize=(12, n_rows*3))

    #fig, ax = plt.subplots(1, 3, figsize=(16, 4))

    for i, species in enumerate(species_list):
        ind_j = i % 3
        ind_i = i // 3
        if n_rows == 1:
            ax_ij = ax[ind_j]
        else:
            ax_ij = ax[ind_i, ind_j]

        sns.boxplot(x='species', y='concentration', hue='model_type', data=data_df[data_df['species'] == species],
                    ax=ax_ij, palette=colorlist_swapped)
        #sns.boxplot(x='species', y='concentration', hue='model_type', data=data_df[data_df['species'] == species],
        #            ax=ax[i], palette=colorlist_swapped)
        ax_ij.legend_.remove()
        ax_ij.set_xlabel('')
        ax_ij.set_ylabel(y_label)
        ax_ij.yaxis.grid(True)

        #ax[i].legend_.remove()
        #ax[i].set_xlabel('')
        #ax[i].set_ylabel(y_label)
        #ax[i].yaxis.grid(True)

        plt.sca(ax_ij)
        #plt.sca(ax[i])
        fig.gca().get_yaxis().set_major_formatter(ScalarFormatter(useOffset=False))
        fig.gca().get_yaxis().set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    plt.tick_params(axis='x',  # changes apply to the x-axis
                    which='both',  # both major and minor ticks are affected
                    bottom='off',  # ticks along the bottom edge are off
                    top='off',  # ticks along the top edge are off
                    labelbottom='off')  # labels along the bottom edge are off

    fig.legend([plt.Line2D((0, 1), (0, 0), color=colorlist_swapped[i], linestyle='-') for i in range(len(model_types))],
               model_types,
               numpoints=1,
               ncol=2,
               loc='center',
               bbox_to_anchor=(0.5, 0.95))

    fig.subplots_adjust(top=0.85, wspace=0.4, hspace=0.5)

    plt.savefig(''.join([file_out_base, '_multiple.pdf']))
    plt.close()

    return fig, ax


def _plot_aggregated_boxplots(median_dic_dKd, median_dic_nodKd, file_out_base, y_lims, y_label, plot_spacings, fig_size,
                              legend_label=None):

    data_df_dKd = pd.DataFrame.from_dict(median_dic_dKd)
    data_df_dKd['model_type'] = ['dKd' for i in range(len(data_df_dKd.index))]

    data_df_nodKd = pd.DataFrame.from_dict(median_dic_nodKd)
    data_df_nodKd['model_type'] = ['nodKd' for i in range(len(data_df_nodKd.index))]

    data_df = pd.concat([data_df_nodKd, data_df_dKd])

    species_list = data_df_nodKd.columns.values[:-1]
    print data_df_nodKd.columns.values
    print species_list
    model_types = ['nodKd', 'dKd']

    df_long = pd.melt(data_df, 'model_type', var_name='species', value_name='concentration')

    colorlist_swapped = [COLOR_LIST[0], COLOR_LIST[1]]  # , COLOR_LIST[2]]

    _plot_single(df_long, model_types, species_list, colorlist_swapped, file_out_base, y_label, legend_label,
                 plot_spacings, fig_size)
    _plot_subplots(df_long, model_types, species_list, colorlist_swapped, file_out_base, y_label)


def _get_single_time_point_data(data_df_nodKd, data_df_dKd, species_list, set_i, time_point_ind):

    data_df = pd.DataFrame()

    n_models_dKd = len(data_df_dKd[species_list[0]].ix[time_point_ind, :].values)
    n_models_nodKd = len(data_df_nodKd[species_list[0]].ix[time_point_ind, :].values)

    for species in species_list:
        data_df[species] = list(data_df_nodKd[species].ix[time_point_ind, :].values) + \
                           list(data_df_dKd[species].ix[time_point_ind, :].values)

    data_df['set_i'] = [set_i for i in range(n_models_nodKd)] + [set_i for i in range(n_models_dKd)]
    data_df['model_type'] = ['nodKd' for i in range(n_models_nodKd)] + ['dKd' for i in range(n_models_dKd)]

    return data_df, n_models_nodKd, n_models_dKd


def _save_dic_as_dataframe(dic, file_out, label):

    data_df_temp = pd.DataFrame.from_dict(dic)
    data_df_temp.to_csv(''.join([file_out, label, '.csv']), sep='\t')


def _load_dataframe_to_dic(file_in, label):

    data_df_temp = pd.read_csv(''.join([file_in, label, '.csv']), sep='\t', index_col=0)
    dic = {}
    for key in data_df_temp.columns:
        dic[key] = data_df_temp[key].values
    return dic


def plot_single_time_point_boxplot(data_df, file_out_base, species_list, y_lims, set_i):
    #seaborn_boxplots()

    n_species = len(species_list)
    model_types = ['nodKd', 'dKd']
    fig, ax = plt.subplots(1, n_species, figsize=(n_species * 4, 4))
    colorlist_swapped = [COLOR_LIST[1], COLOR_LIST[0], COLOR_LIST[2]]
    #

    for species_i, species in enumerate(species_list):
        sns.boxplot(x='set_i', y=species, hue='model_type', data=data_df,  # split=True,
                    palette=colorlist_swapped, ax=ax[species_i])

        ax[species_i].set_xlabel(species)
        ax[species_i].set_ylabel('Concentration (mol/L)')
        if y_lims:
            ax[species_i].set_ylim(y_lims[species_i])
        ax[species_i].legend_.remove()
        ax[species_i].grid()

        plt.sca(ax[species_i])  # set axis current axis active
        plt.tick_params(axis='x',  # changes apply to the x-axis
                        which='both',  # both major and minor ticks are affected
                        bottom='off',  # ticks along the bottom edge are off
                        top='off',  # ticks along the top edge are off
                        labelbottom='off')  # labels along the bottom edge are off

        fig.gca().get_yaxis().set_major_formatter(ScalarFormatter(useOffset=False))

    fig.legend([plt.Line2D((0, 1), (0, 0), color=colorlist_swapped[i], linestyle='-') for i in range(len(model_types))],
               model_types,
               numpoints=1,
               ncol=2,
               loc='center',
               bbox_to_anchor=(0.5, 0.95))

    fig.subplots_adjust(top=0.85, bottom=0.15, left=0.1, right=0.95, wspace=0.4)

    plt.savefig(''.join([file_out_base, str(set_i), '_boxplot.png']))
    plt.close()


'''


def get_species_range(data_df, dic_dKd, dic_nodKd, offset=1):

    for species in data_df.columns.values[:-offset]:
        dic_dKd[species].append(np.max(data_df[data_df['model_type'] == 'dKd'][species]) - np.min(data_df[data_df['model_type'] == 'dKd'][species]))
        dic_nodKd[species].append(np.max(data_df[data_df['model_type'] == 'nodKd'][species]) - np.min(data_df[data_df['model_type'] == 'nodKd'][species]))



def plot_filtered_timecourses(data_df_dic, file_out_base, species_list, set_i, x_lims, y_lims, label):
    """
    Given dataframes with and without dKd, plots data in a figure with 2 subplots.

    :param data_df_nodKd:
    :param data_df_dKd:
    :param file_out_base:
    :param species_list:
    :param set_i:
    :param x_lims:
    :param y_lims:
    :return:
    """

    n_df = len(data_df_dic)
    n_cols = 3
    n_rows = n_df / n_cols if n_df % n_cols == 0 else n_df / n_cols + 1

    fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(18, 4 * n_rows), sharey=True)

    for i, bin in enumerate(data_df_dic.keys()):
        ind_j = i % 3
        ind_i = i // 3
        _subplot_filtered(ax[ind_i, ind_j], data_df_dic[bin], species_list, x_lims, y_lims, bin)

    ax[0, 1].legend([plt.Line2D((0, 1), (0, 0), color=COLOR_LIST[i], linestyle='-') for i in range(len(species_list))],
                    species_list,
                    numpoints=1,
                    ncol=len(species_list),
                    loc='upper center',
                    bbox_to_anchor=(0.5, 1.4))

    fig.subplots_adjust(bottom=0.15, left=0.09, right=0.89, hspace=0.3)
    plt.savefig(''.join([file_out_base, '_', str(set_i), '_', label, '_by_bin.png']))
    plt.close()





def filter_df(data_df, model_list, n_species):
    """
    Given a dataframe with time course data, where there is such data for different models, and given a list of models,
    it creates a new dataframe with time-course data only for the models in model_list.

    :param data_df: dataframe with original time course data
    :param model_list: a list of integers that represent the model number
    :param n_species: how many species are considered (e.g., metabolites, enzymes, fluxes)
    :return: a dataframe with time-course data only for the models in model_list
    """

    column_ranges = [[item * n_species + 1, (item) * n_species + n_species + 1] for item in
                     model_list]  # if not math.isnan(item)]
    data_df_filtered = data_df['t']

    for i, (pos_1, pos_2) in enumerate(column_ranges):
        data_df_filtered = pd.concat([data_df_filtered, data_df.ix[:, pos_1:pos_2]], axis=1)

    return data_df_filtered


def filter_models_by_concentration(data_df, concentration_threshold, time_point_ind, offset, n_species, flag_sign):
    """
    Given a dataframe with time course data, a concentration threshold, time point index, and a flag sign it returns a
    list with the models that satisfy, e.g.:
        concentration at a given time point (represented by its index in dataframe) > concentration threshold,
    if flag_sign='>'.

    :param data_df: dataframe with time-course data
    :param concentration_threshold: a float with threshold for concentration
    :param time_point_ind: the index of the time point in data_df at which to check the concentration value
    :param offset: index of the first column for the concentration of the species of interest
    :param n_species: how many species are represented in data_df
    :param flag_sign: either '>' for condition concentration > concentration_threshold or '<' for concentration < concentration_threshold
    :return:
    """

    model_list = None

    if flag_sign == '>':
        model_list = [(i - offset) / n_species for i in range(offset, len(data_df.columns), n_species)
                      if data_df.ix[time_point_ind, i] > concentration_threshold]

    elif flag_sign == '<':
        model_list = [(i - offset) / n_species for i in range(offset, len(data_df.columns), n_species)
                      if data_df.ix[time_point_ind, i] < concentration_threshold]

    return model_list


def determine_convergence(data_df, t_i, t_f, concentration_threshold, offset, n_species):
    """
    Given a dataframe with time course data, initial and final time points, and a threshold for the difference in
     concentration between those time points, it returns a list with the indices of the models where the concentration
       of a given species at t_f - concentration at t_i is greater or equal to concentration_threshold.

    :param data_df: dataframe with time course data
    :param t_i: initial time point
    :param t_f: final time point
    :param concentration_threshold: concentration above which  the model is considered to not have converged
    :param offset: the first column index of the species concentrations
    :param n_species: how many species are includcd in dataframe
    :return: list with the indexes of the models that do not converge
    """

    ind_i = data_df[data_df['t'] == t_i].index.values[0]
    ind_f = data_df[data_df['t'] == t_f].index.values[0]

    model_list = [(i - offset) / n_species for i in range(offset, len(data_df.columns), n_species)
                  if abs(data_df.ix[ind_f, i] - data_df.ix[ind_i, i]) >= concentration_threshold]

    return model_list


def get_bin_model_list(file_in):
    data_df = pd.read_csv(file_in, sep=',', index_col=0)
    data_dic = data_df.to_dict(orient='list')

    for key in data_df.keys():
        data_dic[key] = [int(value) for value in data_dic[key] if not math.isnan(value)]

    return data_dic

'''


def time_course_analysis(file_in_base, file_out_base, species_list, Keq, subs_list, prod_list, mets_order, legend_label,
                         n_model_ensembles, plot_type, time_point_list, plot_spacings, fig_size, y_lims=None,
                         y_lims_timecourses=None, plot_only=False):
    """

    Analyze model simulations: time-courses for metabolite and enzyme concentrations as well as reactions' fluxes.
    Plots the time-courses and median concentrations/fluxes at t=0 and t=t_final.

    :param file_in_base: path and base name for the input file
    :param file_out_base: path and base name for the output file
    :param species_list: list with species being analyzed: metabolite, enzyme, or reactions names/ids.
    :param Keq: equilibrium constant value
    :param subs_list: list of reaction substrates
    :param prod_list: list of reaction products
    :param mets_order: list with metabolite/enzyme/reaction order
    :param legend_label: labels for the legend
    :param n_model_ensembles: how many model ensembles are being analyzed
    :param plot_type: kind of data to be plotted: metabolite, enzyme, reaction flux
    :param time_point_list: list of time points to be analyzed, tipically t=0, and t=t_final
    :param plot_spacings: space between plots for some boxplots
    :param fig_size: figure size for the box plots
    :param y_lims: tuple or list with limits for y-axis in boxplots
    :param y_lims_timecourses: tuple or list with limits for y-axis in time-course plots
    :param plot_only: boolean specifying whether data is to be analyzed or only plot boxplots
    :return: None
    """

    if plot_type == 'flux':
        median_dic_dKd = {species: [] for species in species_list}
        median_dic_nodKd = {species: [] for species in species_list}

    if plot_type == 'enz' or plot_type == 'mets':
        median_dic_dKd = {}
        median_dic_dKd[time_point_list[0]] = {species: [] for species in species_list}
        median_dic_dKd[time_point_list[1]] = {species: [] for species in species_list}

        median_dic_nodKd = {}
        median_dic_nodKd[time_point_list[0]] = {species: [] for species in species_list}
        median_dic_nodKd[time_point_list[1]] = {species: [] for species in species_list}

    if not plot_only:
        for ensemble_i in range(n_model_ensembles):

            try:
                data_df_nodKd = import_data_timecourses(''.join([file_in_base, '_dKd_', str(ensemble_i), '_', str(Keq), '.csv']),
                                                        Keq, subs_list, prod_list)
                data_df_dKd = import_data_timecourses(''.join([file_in_base, '_all_', str(ensemble_i), '_', str(Keq), '.csv']), Keq,
                                                      subs_list, prod_list)
            except IOError:
                print 'could not find: ' + ''.join([file_in_base, str(ensemble_i), '_', str(Keq), '.csv'])
                continue

            plot_timecourses(data_df_nodKd, data_df_dKd, file_out_base, species_list, ensemble_i, y_lims_timecourses,
                             Keq, mets_order, legend_label)

            if plot_type == 'flux':

                time_point = time_point_list[-1]
                time_point_ind = data_df_dKd[data_df_dKd.ix[:, 0] == time_point].index.values[0]

                data_df, n_models_nodKd, n_models_dKd = _get_single_time_point_data(data_df_nodKd, data_df_dKd,
                                                                                    species_list,
                                                                                    ensemble_i, time_point_ind)

                #plot_single_time_point_boxplot(data_df, ''.join([file_out_base, '_last_t_']),
                #                               species_list, y_lims, ensemble_i)

                get_species_median(data_df, median_dic_dKd, median_dic_nodKd, offset=2)

            if plot_type == 'enz' or plot_type == 'mets':

                for time_point in time_point_list:
                    time_point_ind = data_df_dKd[data_df_dKd.ix[:, 0] == time_point].index.values[0]
                    data_df, n_models_nodKd, n_models_dKd = _get_single_time_point_data(data_df_nodKd, data_df_dKd,
                                                                                        species_list,
                                                                                        ensemble_i, time_point_ind)

                    #plot_single_time_point_boxplot(data_df,
                    #                               ''.join([file_out_base, '_first_t_' if time_point == 0 else '_last_t_']),
                    #                               species_list, y_lims, ensemble_i)

                    get_species_median(data_df, median_dic_dKd[time_point], median_dic_nodKd[time_point], offset=2)

    if plot_type == 'flux':

        file_out_median = ''.join([file_out_base, '_aggregated_boxplot_',  str(Keq), '_', 'last_t'])

        label_median_dKd = 'median_dic_dKd'
        label_median_nodKd = 'median_dic_nodKd'

        if plot_only:
            median_dic_dKd = _load_dataframe_to_dic(file_out_median, label_median_dKd)
            median_dic_nodKd = _load_dataframe_to_dic(file_out_median, label_median_nodKd)

        else:
            for species in median_dic_dKd:
                median_dic_dKd[species] = median_dic_dKd[species]
                median_dic_nodKd[species] = median_dic_nodKd[species]

            _save_dic_as_dataframe(median_dic_dKd, file_out_median, label_median_dKd)
            _save_dic_as_dataframe(median_dic_nodKd, file_out_median, label_median_nodKd)

        y_label = 'Median concentration at steady-state (mol/L)'
        _plot_aggregated_boxplots(median_dic_dKd, median_dic_nodKd, file_out_median, y_lims, y_label, plot_spacings,
                                  fig_size)

    if plot_type == 'enz' or plot_type == 'mets':

        for time_point in time_point_list:

            file_out_median = ''.join([file_out_base, '_aggregated_boxplot_',  str(Keq), '_', 'first_t' if time_point == 0 else 'last_t'])

            label_median_dKd = 'median_dic_dKd'
            label_median_nodKd = 'median_dic_nodKd'

            if plot_only:
                median_dic_dKd[time_point] = _load_dataframe_to_dic(file_out_median, label_median_dKd)
                median_dic_nodKd[time_point] = _load_dataframe_to_dic(file_out_median, label_median_nodKd)

            else:
                _save_dic_as_dataframe(median_dic_dKd[time_point], file_out_median, label_median_dKd)
                _save_dic_as_dataframe(median_dic_nodKd[time_point], file_out_median, label_median_nodKd)

            for species in median_dic_dKd[time_point]:

                max_val = max(max(median_dic_dKd[time_point][species]), max(median_dic_nodKd[time_point][species]))
                median_dic_dKd[time_point][species] = median_dic_dKd[time_point][species] / max_val
                median_dic_nodKd[time_point][species] = median_dic_nodKd[time_point][species] / max_val


            y_label = 'Relative median concentration'
            _plot_aggregated_boxplots(median_dic_dKd[time_point], median_dic_nodKd[time_point], file_out_median, y_lims, y_label,
                                      plot_spacings, fig_size, legend_label)