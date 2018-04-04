
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

from src.kinetics_integration.binning import bin_keqs
from src.kinetics_integration.plots_definitions import entropy_plot_single_defs
from src.kinetics_integration.utils import import_rateconstants

COLOR_LIST = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628']



def merge_data_frames(base_dir, enz_name, dataset, n_ensembles, limit):


    file_in = ''.join([base_dir, enz_name, '_rateconst_fixed_width_3_', dataset, '_', str(0), '_', str(limit), '_pattern_count.csv'])
    data_df_complete = pd.read_csv(file_in, sep=',')

    n_patterns = 0
    for ensemble_i in range(1, n_ensembles):
        #print(ensemble_i)

        file_in = ''.join([base_dir,  enz_name, '_rateconst_fixed_width_3_', dataset, '_', str(ensemble_i), '_', str(limit), '_pattern_count.csv'])
        data_df_temp = pd.read_csv(file_in, sep=',')
        #print data_df_temp
        data_df_complete = data_df_complete.append(data_df_temp, ignore_index=True)
        n_patterns += len(data_df_temp.index)


    data_df_complete.columns = ['patterns', 'count']
    data_df_complete.to_csv(''.join([base_dir, enz_name, '_rateconst_fixed_width_3_', dataset, '_', str(limit), '_all.csv']))
    print dataset
    print 'average ensemble size: ', data_df_complete['count'].sum() / float(n_ensembles), '\n'
    print 'average patterns per ensemble: ', n_patterns / float(n_ensembles), '\n'

    print('n models: ', data_df_complete['count'].sum())
    print('n patterns: ', len(data_df_complete['patterns'].values))
    print('unique patterns: ', len(set(data_df_complete['patterns'].values)))
    print '\n'

    return data_df_complete


def _plot_results(base_dir, enz_name, dataset_list, ind, n_models_list, n_unique_patterns_list, unique_patterns_n_models_ratio, pattern_intersection, boxplot_xlabel):

    entropy_plot_single_defs()
    plt.rcParams['axes.labelpad'] = '10'

    bar_width = 0.25

    fig, ax = plt.subplots(figsize=(6.5, 4.2))

    #rects1 = ax.bar(ind, n_models_list, width=bar_width, color=COLOR_LIST[1])
    #rects2 = ax.bar(np.array(ind) + bar_width, n_unique_patterns_list, width=bar_width, color=COLOR_LIST[2])
    #rects3 = ax.bar(np.array(ind) + 2*bar_width, pattern_intersection, width=bar_width, color=COLOR_LIST[3])
    print np.array(n_unique_patterns_list)
    print np.array(n_models_list) * 100
    print np.array(n_unique_patterns_list)/np.array(n_models_list) * 100
    print(np.array(pattern_intersection)/n_unique_patterns_list[0]) * 100

    rects1 = ax.bar(ind, np.array(n_unique_patterns_list)/np.array(n_models_list) * 100, width=bar_width, color=COLOR_LIST[1])
    rects2 = ax.bar(np.array(ind) + bar_width, np.array(pattern_intersection)/n_unique_patterns_list[0] * 100, width=bar_width, color=COLOR_LIST[2])

    #plt.plot((-bar_width, len(ind)), (pattern_intersection[0], pattern_intersection[0]), linestyle='--', color=COLOR_LIST[0], linewidth=2)

    #ax.set_yscale('log')
    ax.set_xlabel('Data point removed')
    ax.set_xticks(np.array(ind) + 1.*bar_width)
    ax.set_xticklabels(boxplot_xlabel)
    ax.set_xlim(-bar_width, len(ind))
    ax.set_ylim(0, 100)
    ax.grid(True)

    #ax.legend((rects1[0], rects2[0], rects3[0]),
    ax.legend((rects1[0], rects2[0]),
              #('$n$ valid patterns', '$n$ different patterns', '$n$ patterns overlap'),
              ('% different patterns',  ' % overlap'),
               ncol=2,
               loc='center',
               bbox_to_anchor=(0.5, 1.06))


    fig.subplots_adjust(bottom=0.2, top=0.91, left=0.08, right=0.98)
    plt.savefig(''.join([base_dir, '_', enz_name, '_check_patterns_over_all_models.pdf']))
    plt.close()


    fig, ax = plt.subplots()

    plt.bar(ind, unique_patterns_n_models_ratio, color=COLOR_LIST[1], align='center')
    plt.xticks(range(len(dataset_list)), boxplot_xlabel)

    ax.set_xlabel('Data point removed')
    ax.set_ylabel('% of unique patterns over all valid models')
    ax.grid(True)

    plt.tight_layout()
    plt.savefig(''.join([base_dir, '_check_patterns_per_n_models.pdf']))
    plt.close()



def plot_ENO(base_dir):


    enz_name = 'ENO'
    merge_flag = True
    limit = 100

    base_dir = ''.join([base_dir, '/', enz_name, '/', enz_name, '_param_inf2/output/entropy/'])

    boxplot_xlabel = ['None', '$\Delta K_b$', '$K^{2PG}_{m}$',' $k^{f}_{cat}$', '$K_{eq}$' ]

    dataset_list = ['all', 'dKd', 'Km', 'kcat', 'Keq']
    n_ensembles = 100

    if merge_flag:
        data_df_all_complete = merge_data_frames(base_dir, enz_name, 'all', n_ensembles, limit)
    else:
        data_df_all_complete = pd.read_csv(''.join([base_dir, enz_name, '_rateconst_fixed_width_3_', 'all', '_', str(limit), '_all.csv']))

    data_df_list = []
    n_models_list = []
    n_unique_patterns_list = []
    unique_patterns_n_models_ratio = []
    pattern_intersection = []
    for i, dataset in enumerate(dataset_list):

        if merge_flag:
            data_df_list.append(merge_data_frames(base_dir, enz_name, dataset, n_ensembles, limit))
        else:
            data_df_list.append(pd.read_csv(''.join([base_dir, enz_name, '_rateconst_fixed_width_3_', dataset, '_', str(limit), '_all.csv'])))

        n_models_list.append(float(data_df_list[i]['count'].sum()))
        n_unique_patterns_list.append(float(len(set(data_df_list[i]['patterns'].values))))
        unique_patterns_n_models_ratio.append(len(set(data_df_list[i]['patterns'].values)) / float(data_df_list[i]['count'].sum()) * 100)


        pattern_intersection.append(float(len(set.intersection(set(data_df_all_complete['patterns'].values), set(data_df_list[i]['patterns'].values)))))

    print pattern_intersection

    ind = range(len(dataset_list))
    _plot_results(base_dir, enz_name, dataset_list, ind, n_models_list, n_unique_patterns_list, unique_patterns_n_models_ratio, pattern_intersection, boxplot_xlabel)



def plot_GAPD(base_dir):

    enz_name = 'GAPD'
    merge_flag = True
    limit = 100
    base_dir = ''.join([base_dir, '/', enz_name, '/GAPD_param_inf2/output/entropy/'])

    boxplot_xlabel = ['None', '$\Delta K_b$', '$K_{eq}$', '$K_d^{NAD}$', '$K_m^{NAD}$', '$K_m^{G3P}$', '$K_m^{PI}$', '$k_{cat}$']

    dataset_list = ['all',  'dKd', 'Keq', 'Kd', 'Km1', 'Km2', 'Km3', 'kcat']
    n_ensembles = 100

    if merge_flag:
        data_df_all_complete = merge_data_frames(base_dir, enz_name, 'all', n_ensembles, limit)
    else:
        data_df_all_complete = pd.read_csv(''.join([base_dir, enz_name, '_rateconst_fixed_width_3_', 'all', '_', str(limit), '_all.csv']))

    data_df_list = []
    n_models_list = []
    n_unique_patterns_list =[]
    unique_patterns_n_models_ratio = []
    pattern_intersection = []
    for i, dataset in enumerate(dataset_list):

        if merge_flag:
            data_df_list.append(merge_data_frames(base_dir, enz_name, dataset, n_ensembles, limit))
        else:
            data_df_list.append(pd.read_csv(''.join([base_dir, enz_name, '_rateconst_fixed_width_3_', dataset, '_', str(limit), '_all.csv'])))

        n_models_list.append(float(data_df_list[i]['count'].sum()))
        n_unique_patterns_list.append(float(len(set(data_df_list[i]['patterns'].values))))
        unique_patterns_n_models_ratio.append(len(set(data_df_list[i]['patterns'].values)) / float(data_df_list[i]['count'].sum()) * 100)


        pattern_intersection.append(float(len(set.intersection(set(data_df_all_complete['patterns'].values), set(data_df_list[i]['patterns'].values)))))


    print pattern_intersection

    ind = range(len(dataset_list))
    _plot_results(base_dir, enz_name, dataset_list, ind, n_models_list, n_unique_patterns_list, unique_patterns_n_models_ratio, pattern_intersection, boxplot_xlabel)


def plot_TALA2(base_dir):

    enz_name = 'TALA2'
    limit = 100
    merge_flag = True
    base_dir = ''.join([base_dir, '/', enz_name, '/TALA2_param_inf2/output/entropy/'])

    dataset_list = ['all', 'dKd',  'Keq', 'Km1', 'Km2', 'Km3', 'Km4', 'kcat', 'Ki']
    boxplot_xlabel = ['None', '$\Delta K_b$', '$K_{eq}$', '$K_m^{G3P}$', '$K_m^{E4P}$', '$K_m^{F6P}$', '$K_m^{S7P}$', '$k_{cat}$', '$K_i$']


    n_ensembles = 100

    if merge_flag:
        data_df_all_complete  = merge_data_frames(base_dir, enz_name, 'all', n_ensembles, limit)
    else:
        data_df_all_complete = pd.read_csv(''.join([base_dir, enz_name, '_rateconst_fixed_width_3_', 'all', '_', str(limit), '_all.csv']))

    data_df_list = []
    n_models_list = []
    n_unique_patterns_list =[]
    unique_patterns_n_models_ratio = []
    pattern_intersection = []
    for i, dataset in enumerate(dataset_list):

        if merge_flag:
            data_df_list.append(merge_data_frames(base_dir, enz_name, dataset, n_ensembles, limit))
        else:
            data_df_list.append(pd.read_csv(''.join([base_dir, enz_name, '_rateconst_fixed_width_3_', dataset, '_', str(limit), '_all.csv'])))

        n_models_list.append(float(data_df_list[i]['count'].sum()))
        n_unique_patterns_list.append(float(len(set(data_df_list[i]['patterns'].values))))
        unique_patterns_n_models_ratio.append(len(set(data_df_list[i]['patterns'].values)) / float(data_df_list[i]['count'].sum()) * 100)


        pattern_intersection.append(float(len(set.intersection(set(data_df_all_complete['patterns'].values), set(data_df_list[i]['patterns'].values)))))

    print pattern_intersection

    ind = range(len(dataset_list))
    _plot_results(base_dir, enz_name, dataset_list, ind, n_models_list, n_unique_patterns_list, unique_patterns_n_models_ratio, pattern_intersection, boxplot_xlabel)



def main():
    base_dir = '/home/mrama/Desktop/MD/eMASS-MD_complete_data/enzyme_models/'

    plot_ENO(base_dir)
    plot_GAPD(base_dir)
    plot_TALA2(base_dir)




if __name__ == '__main__':
    main()