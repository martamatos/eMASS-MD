import re
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

from src.kinetics_integration.binning import bin_keqs
from src.kinetics_integration.utils import import_rateconstants

COLOR_LIST = sns.color_palette("Set1")


def plot_keq_bin(data_df_keqs, data_df_keq_bins, v_min, v_max, file_out):
    n_bins = len(data_df_keq_bins.columns)
    n_cols = 1
    n_rows = n_bins / n_cols + 1 if n_bins % n_cols > 0 else n_bins / n_cols

    # fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=[3*n_cols, 2*n_rows], sharex=True)
    fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=[3 * n_cols, n_rows], sharex=True)

    bin_i = 0
    for i in range(n_rows):

        for j in range(n_cols):

            if n_rows == 1:
                sns.heatmap(data_df_keqs.ix[data_df_keq_bins.ix[:, bin_i].dropna(), :], vmin=v_min, vmax=v_max,
                            ax=ax[j])
                plt.sca(ax[j])  # set axis current axis active
            elif n_cols == 1:
                sns.heatmap(data_df_keqs.ix[data_df_keq_bins.ix[:, bin_i].dropna(), :], vmin=v_min, vmax=v_max,
                            ax=ax[i])
                plt.sca(ax[i])  # set axis current axis active
            else:
                sns.heatmap(data_df_keqs.ix[data_df_keq_bins.ix[:, bin_i].dropna(), :], vmin=v_min, vmax=v_max,
                            ax=ax[i, j])
                plt.sca(ax[i, j])  # set axis current axis active

            plt.title(data_df_keq_bins.columns.values[bin_i])
            bin_i += 1
            if bin_i == n_bins:
                break

    plt.tight_layout()
    plt.savefig(''.join([file_out, '_keq_bin.png']))
    plt.close()


def _get_common_stats(pattern_df, f_out):
    mean = pattern_df.ix[:, 0].mean()
    std = pattern_df.ix[:, 0].std(ddof=1)

    f_out.write(''.join(['n models:\t', str(pattern_df.ix[:, 0].sum()), '\n']))
    f_out.write(''.join(['n bins:\t', str(len(pattern_df.index)), '\n']))
    f_out.write(''.join(['models per bin - min:\t%.2f' % pattern_df.ix[:, 0].min(), '\n']))
    f_out.write(''.join(['models per bin - max:\t%.2f' % pattern_df.ix[:, 0].max(), '\n']))
    f_out.write(''.join(['models per bin - average:\t%.2f' % mean, '\n']))
    f_out.write(''.join(['models per bin - std:\t%.2f' % std, '\n']))
    f_out.write(''.join(['models per bin - std in percent:\t%.2f' % ((std / mean) * 100), '\n']))

    conf_int = stats.norm.interval(0.95, loc=mean, scale=std)
    f_out.write(''.join(['95% confidence interval:\t', str(conf_int), '\n']))


def _plot_statistics(plot_data, pattern_df_dic, file_out):

    fig, ax = plt.subplots(1, figsize=[6, 4])

    sns.boxplot(data=plot_data, palette=COLOR_LIST)
    sns.stripplot(data=plot_data, jitter=True, color='black')

    ax.set_xticklabels(pattern_df_dic.keys(), rotation=90)

    ax.legend([plt.Line2D((0, 1), (0, 0), color=COLOR_LIST[1], linestyle='-'),
               plt.Line2D((0, 1), (0, 0), color=COLOR_LIST[0], linestyle='-')],
              ['nodKd', 'dKd'],
              ncol=2,
              numpoints=1,
              loc='upper center',
              bbox_to_anchor=(0.5, 1.2))

    fig.subplots_adjust(top=0.85, bottom=0.2)

    plt.savefig(''.join([file_out, '.png']))
    plt.close()


def gather_statistics(pattern_df_dic, model_type_ref, file_out):

    plot_data = []
    bins_dic = OrderedDict((key, []) for key in pattern_df_dic.keys())

    with open(''.join([file_out, '.csv']), 'w+') as f_out:

        for key in pattern_df_dic:
            f_out.write(''.join(['\n\nkey: ', key, '\n']))
            _get_common_stats(pattern_df_dic[key], f_out)
            bins_dic[key] = set(pattern_df_dic[key].index.values)

            f_out.write(''.join([key, ' bins:\t', str(bins_dic[key]), '\n']))

            if key != model_type_ref:
                n_bins_common = len(bins_dic[key].intersection(bins_dic[model_type_ref]))
                n_bins_unique = len(bins_dic[key].difference(bins_dic[model_type_ref]))
                n_bins = n_bins_common + n_bins_unique

                f_out.write(''.join(
                    ['n common bins with ', model_type_ref, ' in %:\t', str(n_bins_common / float(n_bins)), '\n']))
                f_out.write(''.join(['n unique bins in %:\t', str(n_bins_unique / float(n_bins)), '\n']))

            plot_data.append(pattern_df_dic[key].ix[:, 0].values)

    #_plot_statistics(plot_data, pattern_df_dic, file_out)
    return bins_dic


def _parse_model_type_specific_info(f_in, model_type, model_type_ref, n_models, n_bins, n_models_per_bin, n_common_bins,
                                    n_unique_bins):
    line = f_in.readline()
    n_models_temp = float(re.findall('n models:\t(\d+)', line)[0])
    n_models[model_type].append(n_models_temp)

    line = f_in.readline()
    n_bins_temp = float(re.findall('n bins:\t(\d+)', line)[0])
    n_bins[model_type].append(n_bins_temp)

    f_in.readline()
    f_in.readline()
    line = f_in.readline()
    n_models_per_bin_temp = float(re.findall('models per bin - average:\t(\d+.\d+)', line)[0])
    n_models_per_bin[model_type].append(n_models_per_bin_temp)

    if model_type != model_type_ref:
        f_in.readline()
        f_in.readline()
        f_in.readline()
        f_in.readline()
        line = f_in.readline()
        n_common_bins_temp = float(
            re.findall(''.join(['n common bins with ', model_type_ref, ' in %:\t(\d+.\d+)']), line)[0])
        n_common_bins[model_type].append(n_common_bins_temp)

        line = f_in.readline()
        n_unique_bins_temp = float(re.findall('n unique bins in %:\t(\d+.\d+)', line)[0])
        n_unique_bins[model_type].append(n_unique_bins_temp)


def _calculate_sample_means_stats(data_list, f_out, list_name):
    min = np.array(data_list).min()
    max = np.array(data_list).max()
    mean = np.array(data_list).mean()
    std = np.array(data_list).std() / np.sqrt(len(data_list))
    conf_int = stats.norm.interval(0.95, loc=mean, scale=std)

    f_out.write(''.join([list_name, '\n']))
    f_out.write(''.join(['mean: %.2f' % mean, '\n']))
    f_out.write(''.join(['std: %.2f' % std, '\n']))
    f_out.write(''.join(['std in percentage: %.2f' % ((std / mean) * 100), '\n']))
    f_out.write(''.join(['min: %.2f' % min, '\n']))
    f_out.write(''.join(['max: %.2f' % max, '\n']))
    f_out.write(''.join(['confidence interval', str(conf_int), '\n\n']))

    return conf_int


def _plot_aggregated_statistics(n_models_data, n_bins_data, n_models_per_bin_data, n_common_bins_dic_data,
                                n_unique_bins_dic_data, model_type_list, file_out):

    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=[8, 10])

    sns.boxplot(data=n_models_data, ax=ax[0, 0])
    sns.boxplot(data=n_bins_data, ax=ax[0, 1])
    sns.boxplot(data=n_models_per_bin_data, ax=ax[1, 0])

    sns.boxplot(data=n_common_bins_dic_data, ax=ax[2, 0])
    sns.boxplot(data=n_unique_bins_dic_data, ax=ax[2, 1])

    # means = [(lb+ub) / 2. for lb, ub in confidence_intervals]
    # error = [(ub-lb) / 2. for lb, ub in confidence_intervals]
    ax[0, 0].set_title('n models')
    ax[0, 1].set_title('n bins')
    ax[1, 0].set_title('n models per bin')
    ax[2, 0].set_title('n common bins with ref')
    ax[2, 1].set_title('n unique bins')

    # ax[0, 0].set_ylim(30, 90)
    # ax[0, 1].set_ylim(5, 65)
    # ax[1, 0].set_ylim(1, 7)
    ax[2, 0].set_ylim(0, 1)
    ax[2, 1].set_ylim(0, 1)

    ax[0, 0].set_xticklabels(model_type_list, rotation=90)
    ax[0, 1].set_xticklabels(model_type_list, rotation=90)
    ax[1, 0].set_xticklabels(model_type_list, rotation=90)
    ax[2, 0].set_xticklabels(model_type_list, rotation=90)
    ax[2, 1].set_xticklabels(model_type_list, rotation=90)

    plt.tight_layout()
    plt.savefig(''.join([file_out, '.pdf']))
    plt.close()


def gather_statistics_sample_means(file_in_base, model_type_list, model_type_ref, n_model_sets, bin_width, file_out):

    n_models_dic = OrderedDict((key, []) for key in model_type_list)
    n_bins_dic = OrderedDict((key, []) for key in model_type_list)
    n_models_per_bin_dic = OrderedDict((key, []) for key in model_type_list)
    n_common_bins_dic = OrderedDict((key, []) for key in model_type_list)
    n_unique_bins_dic = OrderedDict((key, []) for key in model_type_list)

    for model_i in range(1, n_model_sets + 1):
        file_in = ''.join([file_in_base, str(bin_width), '_', str(model_i), '_stats.csv'])

        model_type_i = 0
        with open(file_in, 'r') as f_in:
            line = f_in.readline()
            while line:

                if re.match(''.join(['key: ', model_type_list[model_type_i], '\n']), line):
                    _parse_model_type_specific_info(f_in, model_type_list[model_type_i], model_type_ref, n_models_dic,
                                                    n_bins_dic, n_models_per_bin_dic, n_common_bins_dic,
                                                    n_unique_bins_dic)
                    model_type_i += 1
                    if model_type_i == len(model_type_list):
                        break

                line = f_in.readline()

    confidence_intervals = []
    print n_models_dic
    print n_bins_dic
    print n_models_per_bin_dic
    print n_common_bins_dic
    print n_unique_bins_dic

    with open(''.join([file_out, '.dat']), 'w+') as f_out:
        for model_type in model_type_list:
            conf_int = _calculate_sample_means_stats(n_models_dic[model_type], f_out,
                                                     ''.join(['n models per ensemble - ', model_type]))
            confidence_intervals.append(conf_int)

            f_out.write('\n')
            conf_int = _calculate_sample_means_stats(n_bins_dic[model_type], f_out,
                                                     ''.join(['n bins per ensemble - ', model_type]))
            confidence_intervals.append(conf_int)

            f_out.write('\n')
            conf_int = _calculate_sample_means_stats(n_models_per_bin_dic[model_type], f_out,
                                                     ''.join(['n models per bin per ensemble - ', model_type]))
            confidence_intervals.append(conf_int)

            if model_type != model_type_ref:
                f_out.write('\n')
                conf_int = _calculate_sample_means_stats(n_common_bins_dic[model_type], f_out, ''.join(
                    ['n common bins with ', model_type_ref, ' - ', model_type]))
                confidence_intervals.append(conf_int)

                f_out.write('\n')
                conf_int = _calculate_sample_means_stats(n_unique_bins_dic[model_type], f_out,
                                                         ''.join(['n unique bins - ', model_type]))
                confidence_intervals.append(conf_int)

    n_models_data = [n_models_dic[model_type] for model_type in model_type_list]
    n_bins_data = [n_bins_dic[model_type] for model_type in model_type_list]
    n_models_per_bin_data = [n_models_per_bin_dic[model_type] for model_type in model_type_list]
    n_common_bins_dic_data = [n_common_bins_dic[model_type] for model_type in model_type_list]
    n_unique_bins_dic_data = [n_unique_bins_dic[model_type] for model_type in model_type_list]

    _plot_aggregated_statistics(n_models_data, n_bins_data, n_models_per_bin_data, n_common_bins_dic_data,
                                n_unique_bins_dic_data, model_type_list, file_out)


def _write_extra_stats(file_out, bins_dic, n_models_count):

    with open(file_out, 'w') as f_out:

        f_out.write(''.join(['total nodKd models: ', str(n_models_count['nodKd']), '\n']))
        f_out.write(''.join(['total dKd models:\n', str(n_models_count['dKd']), '\n\n']))

        all_nodKd_bins = set.union(*bins_dic['nodKd'])
        f_out.write('union of all nodKd bins\n')
        f_out.write(''.join(['n bins: ', str(len(all_nodKd_bins)), '\n']))
        f_out.write(''.join(['bins:\n', str(all_nodKd_bins), '\n\n']))

        common_nodKd_bins = set.intersection(*bins_dic['nodKd'])
        f_out.write('intersection of all nodKd bins\n')
        f_out.write(''.join(['n bins: ', str(len(common_nodKd_bins)), '\n']))
        f_out.write(''.join(['bins: ', str(common_nodKd_bins), '\n\n']))


        all_dKd_bins = set.union(*bins_dic['dKd'])
        f_out.write('union of all dKd bins\n')
        f_out.write(''.join(['n bins: ', str(len(all_dKd_bins)), '\n']))
        f_out.write(''.join(['bins:\n', str(all_dKd_bins), '\n\n']))

        common_dKd_bins = set.intersection(*bins_dic['dKd'])
        f_out.write('intersection of all dKd bins\n')
        f_out.write(''.join(['n bins: ', str(len(common_dKd_bins)), '\n']))
        f_out.write(''.join(['bins: ', str(common_dKd_bins), '\n\n']))


        intersection_nodKd_dKd = set.intersection(all_nodKd_bins, all_dKd_bins)
        f_out.write('intersection of all nodKd bins with all dKd bins\n')
        f_out.write(''.join(['n bins: ', str(len(intersection_nodKd_dKd)), '\n']))
        f_out.write(''.join(['bins: ', str(intersection_nodKd_dKd), '\n\n']))


        all_nodKd_minus_dKd = all_nodKd_bins.difference(all_dKd_bins)
        f_out.write('all nodKd bins minus all dKd bins\n')
        f_out.write(''.join(['n bins: ', str(len(all_nodKd_minus_dKd)), '\n']))
        f_out.write(''.join(['bins: ', str(all_nodKd_minus_dKd), '\n\n']))

        all_dKd_minus_nodkd = all_dKd_bins.difference(all_nodKd_bins)
        f_out.write('all dKd bins minus all nodKd bins\n')
        f_out.write(''.join(['n bins: ', str(len(all_dKd_minus_nodkd)), '\n']))
        f_out.write(''.join(['bins: ', str(all_dKd_minus_nodkd), '\n\n']))



def bin_enolase(ssd_threshold, scale_data, filter_description=None, model_list_nodKd=None, model_list_dKd=None):
    column_order = [0, 2, 1]
    convert_to_ratios = True
    v_min = -9
    v_max = 9
    ylim = 80


    base_dir = '/home/mrama/Dropbox/Kinetics/Enzymes_new/ENO_param_influence/test_kcat/'
    file_in_base = ''.join([base_dir, 'output/treated_data/rateconst_ENO_'])
    file_out_base = ''.join([base_dir, 'output/binning'])

    model_types = ['nodKd', 'dKd']  # , 'KeqdKd', 'Km1dKd', 'Km2dKd', 'Km3dKd', 'Km4dKd', 'kcat1dKd', 'kcat2dKd', 'kcat1kcat2dKd', 'KidKd']
    # model_types = ['Keq', 'dKd',  'Km1', 'Km2', 'Km3', 'Km4', 'kcat1', 'kcat2', 'kcat1kcat2', 'Ki']
    model_type_ref = 'nodKd'
    n_model_ensembles = 100
    Keq = 5.19

    bin_width_list = [5]
    bins_dic = OrderedDict((key, []) for key in model_types)
    n_models_count = OrderedDict((key, 0) for key in model_types)


    for bin_width in bin_width_list:

        for ensemble_i in range(1, n_model_ensembles + 1):
            pattern_df_dic = OrderedDict()
            for model_type in model_types:
                print model_type

                file_in = ''.join([file_in_base, model_type, '_', str(ensemble_i), '.csv'])
                print file_in

                data_df, fitness = import_rateconstants(file_in, ssd_threshold, Keq=Keq, column_order=column_order,
                                                        convert_to_ratios=convert_to_ratios, filter=True)

                file_out = ''.join(
                    [file_out_base, '/ENO_fixed_width_', str(bin_width), '_', model_type, '_', str(ensemble_i),
                     '_onecolumn'])
                pattern_count_df, pattern_df_by_model = bin_keqs(data_df, bin_width=bin_width)

                pattern_df_by_model.to_csv(
                    ''.join([file_out_base, '/ENO_fixed_width_', str(bin_width), '_', model_type, '_', str(ensemble_i)]))

                pattern_df_dic[model_type] = pattern_count_df

                #plot_keq_bin(data_df, pattern_df_by_model, v_min, v_max, file_out)

                file_out = ''.join([file_out_base, '/GAPD_bin_width_', str(bin_width), '_', str(ensemble_i), '_stats'])
                bins_dic_temp = gather_statistics(pattern_df_dic, model_type_ref, file_out)

                bins_dic[model_type].append(bins_dic_temp[model_type])
                n_models_count[model_type] += len(data_df.index)

        file_out = ''.join([file_out_base, '/ENO_bin_width_', str(bin_width), '_extra_stats'])
        _write_extra_stats(file_out, bins_dic, n_models_count)

        file_out = ''.join([file_out_base, '/ENO_bin_width_', str(bin_width), '_stats_sumup'])
        file_in = ''.join([file_out_base,  '/ENO_bin_width_'])
        gather_statistics_sample_means(file_in, model_types, model_type_ref, n_model_ensembles, bin_width, file_out)


def bin_gapd(ssd_threshold, scale_data, filter_description=None, model_list_nodKd=None, model_list_dKd=None):
    column_order = [0, 2, 4, 5, 3, 1]
    convert_to_ratios = True
    v_min = -9
    v_max = 9
    ylim = 80
    Keq=0.408
    #Keq=None

    base_dir = '/home/mrama/Dropbox/Kinetics/Enzymes_new/GAPD_param_influence/param_influence/'
    file_in_base = ''.join([base_dir, 'output/treated_data/rateconst_GAPD_'])
    file_out_base = ''.join([base_dir, 'output/binning'])

    model_types = ['nodKd', 'dKd']  # , 'KeqdKd', 'Km1dKd', 'Km2dKd', 'Km3dKd', 'Km4dKd', 'kcat1dKd', 'kcat2dKd', 'kcat1kcat2dKd', 'KidKd']
    # model_types = ['Keq', 'dKd',  'Km1', 'Km2', 'Km3', 'Km4', 'kcat1', 'kcat2', 'kcat1kcat2', 'Ki']
    model_type_ref = 'nodKd'
    n_model_ensembles = 100
    limit = 57

    bin_width_list = [5]
    bins_dic = OrderedDict((key, []) for key in model_types)
    n_models_count = OrderedDict((key, 0) for key in model_types)

    for bin_width in bin_width_list:

        for ensemble_i in range(1, n_model_ensembles+1):
            pattern_df_dic = OrderedDict()
            for model_type in model_types:
                print model_type

                file_in = ''.join([file_in_base, model_type, '_', str(ensemble_i), '.csv'])
                print file_in

                data_df, fitness = import_rateconstants(file_in, ssd_threshold, Keq=Keq, column_order=column_order,
                                                        convert_to_ratios=convert_to_ratios, filter=True, limit=limit)

                file_out = ''.join(
                    [file_out_base, '/GAPD_fixed_width_', str(bin_width), '_', model_type, '_', str(ensemble_i),
                     '_onecolumn'])
                pattern_count_df, pattern_df_by_model = bin_keqs(data_df, bin_width=bin_width)

                pattern_df_by_model.to_csv(
                    ''.join([file_out_base, '/GAPD_fixed_width_', str(bin_width), '_', model_type, '_', str(ensemble_i)]))

                pattern_df_dic[model_type] = pattern_count_df

                #plot_keq_bin(data_df, pattern_df_by_model, v_min, v_max, file_out)

                file_out = ''.join([file_out_base, '/GAPD_bin_width_', str(bin_width), '_', str(ensemble_i), '_stats'])
                bins_dic_temp = gather_statistics(pattern_df_dic, model_type_ref, file_out)

                bins_dic[model_type].append(bins_dic_temp[model_type])
                n_models_count[model_type] += len(data_df.index)

        file_out = ''.join([file_out_base, '/GAPD_bin_width_', str(bin_width), '_extra_stats'])
        _write_extra_stats(file_out, bins_dic, n_models_count)

        file_in = ''.join([file_out_base,  '/GAPD_bin_width_'])
        file_out = ''.join([file_out_base, '/GAPD_bin_width_', str(bin_width), '_stats_sumup'])
        gather_statistics_sample_means(file_in, model_types, model_type_ref, n_model_ensembles, bin_width, file_out)


def bin_talb(ssd_threshold, scale_data, filter_description=None, model_list_nodKd=None, model_list_dKd=None):
    column_order = [0, 3, 5, 2, 4, 1]
    convert_to_ratios = True
    v_min = -9
    v_max = 9
    ylim = 30
    Keq = 1/1.19

    base_dir = '/home/mrama/Dropbox/Kinetics/Enzymes_new/TALA2_param_influence/'
    file_in_base = ''.join([base_dir, 'results/treated_data/rateconst_TALA2_'])
    file_out_base = ''.join([base_dir, 'results/binning'])

    model_types = ['nodKd', 'dKd']  # , 'KeqdKd', 'Km1dKd', 'Km2dKd', 'Km3dKd', 'Km4dKd', 'kcat1dKd', 'kcat2dKd', 'kcat1kcat2dKd', 'KidKd']
    # model_types = ['Keq', 'dKd',  'Km1', 'Km2', 'Km3', 'Km4', 'kcat1', 'kcat2', 'kcat1kcat2', 'Ki']
    model_type_ref = 'nodKd'
    n_model_ensembles = 99
    limit = 39

    bin_width_list = [5]
    bins_dic = OrderedDict((key, []) for key in model_types)
    n_models_count = OrderedDict((key, 0) for key in model_types)

    for bin_width in bin_width_list:

        for ensemble_i in range(1, n_model_ensembles+1):
            pattern_df_dic = OrderedDict()
            for model_type in model_types:
                print model_type

                file_in = ''.join([file_in_base, model_type, '_', str(ensemble_i), '.csv'])
                print file_in

                data_df, fitness = import_rateconstants(file_in, ssd_threshold, Keq=Keq, column_order=column_order,
                                                        convert_to_ratios=convert_to_ratios, filter=True, limit=limit)

                file_out = ''.join(
                    [file_out_base, '/TALA2_fixed_width_', str(bin_width), '_', model_type, '_', str(ensemble_i),
                     '_onecolumn'])
                pattern_count_df, pattern_df_by_model = bin_keqs(data_df, bin_width=bin_width)

                pattern_df_by_model.to_csv(
                    ''.join([file_out_base, '/TALA2_fixed_width_', str(bin_width), '_', model_type, '_', str(ensemble_i)]))

                pattern_df_dic[model_type] = pattern_count_df

                #plot_keq_bin(data_df, pattern_df_by_model, v_min, v_max, file_out)

                file_out = ''.join([file_out_base, '/TALA2_bin_width_', str(bin_width), '_', str(ensemble_i), '_stats'])
                bins_dic_temp = gather_statistics(pattern_df_dic, model_type_ref, file_out)

                bins_dic[model_type].append(bins_dic_temp[model_type])
                n_models_count[model_type] += len(data_df.index)

        file_out = ''.join([file_out_base, '/TALA2_bin_width_', str(bin_width), '_extra_stats'])
        _write_extra_stats(file_out, bins_dic, n_models_count)

        file_in = ''.join([file_out_base,  '/TALA2_bin_width_'])
        file_out = ''.join([file_out_base, '/TALA2_bin_width_', str(bin_width), '_stats_sumup'])
        gather_statistics_sample_means(file_in, model_types, model_type_ref, n_model_ensembles, bin_width, file_out)


if __name__ == '__main__':
    ssd_threshold = 1
    scale_data = False

    #bin_talb(ssd_threshold, scale_data)
    #bin_enolase(ssd_threshold, scale_data)
    bin_gapd(ssd_threshold, scale_data)