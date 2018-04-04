
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from src.kinetics_integration.binning import bin_keqs
from src.kinetics_integration.plots_definitions import entropy_plot_single_defs
from src.kinetics_integration.utils import import_rateconstants

COLOR_LIST = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628']


def plot_param_prediction(data_df, file_out, boxplot_xlabel, ssd_threshold, filter=False):

    entropy_plot_single_defs()

    if filter:
        data_df = data_df[data_df['ssd'] < ssd_threshold]
        fitness = data_df['ssd'].values
    del data_df['ssd']
    data_df = data_df.ix[1:, :]

    data_df.columns = boxplot_xlabel
    print(len(data_df.index))

    #fig, ax = plt.subplots(1, 1, figsize=(6.5, 4.2))
    fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))
    #fig, ax = plt.subplots(1, 1, figsize=(10, 4.2))
    sns.boxplot(data=data_df, orient='v', palette=COLOR_LIST, ax=ax)
    plt.xticks(rotation=90)


    ax.set_ylabel('Prediction error (%)')
    #ax.set_ylim(-1, 200)  # FBP2
    ax.set_ylim(-1, 50)  # others

    #fig.subplots_adjust(bottom=0.27, top=0.9, left=0.15, right=0.96)
    plt.tight_layout()
    plt.savefig(''.join([file_out, '.pdf']))
    plt.close()

    """
    fig, ax = plt.subplots(3, 1, figsize=(8, 12), sharex=True, sharey=True)
    data_df['Name'] =data_df.index
    parallel_coordinates(data_df.ix[:5,:], 'Name', ax=ax[0], color=COLOR_LIST)
    parallel_coordinates(data_df.ix[6:10,:], 'Name', ax=ax[1], color=COLOR_LIST)
    parallel_coordinates(data_df.ix[11:15,:], 'Name', ax=ax[2], color=COLOR_LIST)
    for i in range(3):
        ax[i].set_yscale('log')
        ax[i].legend_.remove()
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(''.join([file_out, '_pc.pdf']))
    plt.close()"""

    return 0


def check_pattern_overlap(patterns_list, file_out):

    with open(file_out, 'w+') as f_out:

        pattern_intersection = set.intersection(*patterns_list)
        f_out.write(''.join(['common patterns: ', str(len(pattern_intersection)), '\n\n']))

        for i, pattern in enumerate(patterns_list):
            f_out.write(''.join(['n patterns on ensemble ', str(i), ': ', str(len(pattern)), '\n']))



def merge_data_frames(base_dir, dataset, n_ensembles, fit_flag):

    file_in = ''.join([base_dir, 'predicted_params_error_distribution_', dataset, '_', str(0), '_', fit_flag, '.csv'])
    data_df_complete = pd.read_csv(file_in, sep=',')

    for ensemble_i in range(1, n_ensembles):

        file_in = ''.join([base_dir, 'predicted_params_error_distribution_', dataset, '_', str(ensemble_i), '_', fit_flag, '.csv'])
        data_df_temp = pd.read_csv(file_in, sep=',')
        data_df_complete = data_df_complete.append(data_df_temp.ix[1:, :], ignore_index=True)

    data_df_complete.to_csv(''.join([base_dir, 'predicted_params_error_distribution_', dataset, '_', fit_flag, '.csv']))

    return data_df_complete



def plot_ENO(base_dir):

    enz_name = 'ENO'
    merge_flag = True
    ssd_threshold = 0.1

    base_dir = ''.join([base_dir, '/', enz_name, '/ENO_param_inf2/output/treated_data/'])

    boxplot_xlabel = ['$K^{2PG}_{m}$',' $k^{f}_{cat}$', '$K_{eq}$', '$\Delta K_b$']

    custom_rate_const_labels = ['$\overleftarrow{k_{1}}$', '$\overrightarrow{k_{1}}$',
                                '$\overleftarrow{k_{2}}$', '$\overrightarrow{k_{2}}$',
                                '$\overleftarrow{k_{3}}$', '$\overrightarrow{k_{3}}$']

    dataset_list = ['all', 'Km', 'kcat',  'dKd', 'Keq']
    n_ensembles = 100
    fit_flag = 'log_ssd'


    patterns_list = []
    for dataset in dataset_list:
        print(dataset)

        if merge_flag:
            data_df_complete = merge_data_frames(base_dir, dataset, n_ensembles, fit_flag)
        else:
            data_df_complete = pd.read_csv(''.join([base_dir, 'predicted_params_error_distribution_', dataset, '_', fit_flag, '.csv']))


        for ensemble_i in range(n_ensembles):

            file_in = ''.join([base_dir, 'rateconst_', enz_name, '_', dataset, '_', str(ensemble_i), '.csv'])
            column_order = [0, 1, 4, 5, 2, 3]

            bin_width = 3
            data_df, fitness = import_rateconstants(file_in, ssd_threshold, column_order, convert_to_ratios=False)

            pattern_count_df, pattern_df_by_model = bin_keqs(data_df, bin_width=bin_width, lb=-6, ub=9)

            patterns_list.append(set(pattern_count_df.index))

        file_out = ''.join([base_dir, 'pattern_overlap_', dataset, '.txt'])
        check_pattern_overlap(patterns_list, file_out)

        file_out = ''.join([base_dir, 'predicted_params_error_distribution_', dataset, '_abs_ssd'])
        plot_param_prediction(data_df_complete, file_out, boxplot_xlabel, ssd_threshold, filter=True)




def plot_GAPD(base_dir):

    enz_name = 'GAPD'
    merge_flag = True
    fit_flag = 'log_ssd'
    ssd_threshold = 0.1
    base_dir = ''.join([base_dir, '/', enz_name, '/GAPD_param_inf2/output/treated_data/'])

    boxplot_xlabel = ['$K^{NAD}_{m}$','$K^{G3P}_{m}$','$K^{PI}_{m}$',' $k^{f}_{cat}$', '$K_{eq}$', '$K_{d}$', '$\Delta K_b$']
    custom_rate_const_labels = ['$\overleftarrow{k_{1}}$', '$\overrightarrow{k_{1}}$',
                                '$\overleftarrow{k_{2}}$', '$\overrightarrow{k_{2}}$',
                                '$\overleftarrow{k_{3}}$', '$\overrightarrow{k_{3}}$',
                                '$\overleftarrow{k_{4}}$', '$\overrightarrow{k_{4}}$',
                                '$\overleftarrow{k_{5}}$', '$\overrightarrow{k_{5}}$',
                                '$\overleftarrow{k_{6}}$', '$\overrightarrow{k_{6}}$']

    dataset_list = ['all', 'Keq', 'Km1', 'Km2', 'Km3', 'kcat', 'Kd', 'dKd']
    n_ensembles = 100

    patterns_list = []
    for dataset in dataset_list:
        print(dataset)

        if merge_flag:
            data_df_complete = merge_data_frames(base_dir, dataset, n_ensembles, fit_flag)
        else:
            data_df_complete = pd.read_csv(''.join([base_dir, 'predicted_params_error_distribution_', dataset, '_', fit_flag, '.csv']))

        for ensemble_i in range(1, n_ensembles):

            file_in = ''.join([base_dir, 'rateconst_', enz_name, '_', dataset, '_', str(ensemble_i), '.csv'])
            column_order = [0, 1, 4, 5, 8, 9, 10, 11, 6, 7, 2, 3]

            bin_width = 3
            data_df, fitness = import_rateconstants(file_in, ssd_threshold, column_order, convert_to_ratios=False)

            pattern_count_df, pattern_df_by_model = bin_keqs(data_df, bin_width=bin_width, lb=-6, ub=9)

            patterns_list.append(set(pattern_count_df.index))

        file_out = ''.join([base_dir, 'pattern_overlap_', dataset, '.txt'])
        check_pattern_overlap(patterns_list, file_out)

        file_out = ''.join([base_dir, 'predicted_params_error_distribution_', dataset, '_', fit_flag])
        plot_param_prediction(data_df_complete, file_out, boxplot_xlabel, ssd_threshold, filter=True)



def plot_TALA2(base_dir):

    enz_name = 'TALA2'
    merge_flag = True
    fit_flag = 'log_ssd'
    ssd_threshold = 0.1
    base_dir = ''.join([base_dir, '/', enz_name, '/TALA2_param_inf2/output/treated_data/'])

    boxplot_xlabel = ['$K^{G3P}_{m}$','$K^{E4P}_{m}$','$K^{F6P}_{m}$','$K^{S7P}_{m}$',' $k^{f}_{cat1}$', '$k^{f}_{cat2}$',
                      '$K_{eq}$', '$K_{i1}$', '$K_{i2}$', '$\Delta K_b$']

    custom_rate_const_labels = ['$\overleftarrow{k_{1}}$', '$\overrightarrow{k_{1}}$',
                                '$\overleftarrow{k_{2}}$', '$\overrightarrow{k_{2}}$',
                                '$\overleftarrow{k_{3}}$', '$\overrightarrow{k_{3}}$',
                                '$\overleftarrow{k_{4}}$', '$\overrightarrow{k_{4}}$',
                                '$\overleftarrow{k_{5}}$', '$\overrightarrow{k_{5}}$',
                                '$\overleftarrow{k_{6}}$', '$\overrightarrow{k_{6}}$']

    dataset_list = ['all', 'dKd', 'Km1', 'Km2', 'Km3', 'Km4', 'kcat', 'Ki', 'Keq']
    n_ensembles = 100

    patterns_list = []
    for dataset in dataset_list:
        print(dataset)

        if merge_flag:
            data_df_complete = merge_data_frames(base_dir, dataset, n_ensembles, fit_flag)
        else:
            data_df_complete = pd.read_csv(''.join([base_dir, 'predicted_params_error_distribution_', dataset, '_', fit_flag, '.csv']))

        for ensemble_i in range(n_ensembles):

            file_in = ''.join([base_dir, 'rateconst_', enz_name, '_', dataset, '_', str(ensemble_i), '.csv'])
            column_order = [0, 1, 6, 7, 10, 11, 4, 5, 8, 9, 2, 3]

            bin_width = 3
            data_df, fitness = import_rateconstants(file_in, ssd_threshold, column_order, convert_to_ratios=False)

            pattern_count_df, pattern_df_by_model = bin_keqs(data_df, bin_width=bin_width, lb=-6, ub=9)

            patterns_list.append(set(pattern_count_df.index))

        file_out = ''.join([base_dir, 'pattern_overlap_', dataset, '.txt'])
        check_pattern_overlap(patterns_list, file_out)

        file_out = ''.join([base_dir, 'predicted_params_error_distribution_', dataset, '_', fit_flag])
        plot_param_prediction(data_df_complete, file_out, boxplot_xlabel, ssd_threshold, filter=True)


def main():
    base_dir = '/home/mrama/Desktop/MD/eMASS-MD_complete_data/enzyme_models/'
    #base_dir = '/home/mrama/Desktop/MD/eMASS-MD_complete_data/enzyme_models/ENO/ENO_param_scan/output/treated_data/predicted_params_error_distribution_0.1_1_abs_ssd.csv'

    #plot_ENO(base_dir)
    #plot_GAPD(base_dir)
    plot_TALA2(base_dir)




if __name__ == '__main__':
    main()