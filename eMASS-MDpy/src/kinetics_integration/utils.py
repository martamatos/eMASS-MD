import random

import numpy as np
import pandas as pd
from scipy import stats


def set_df_column_labels(data_df, column_order=None, convert_to_ratios=False):
    """
    Sets the column names for data_df, a pandas dataframe.
    If column_order is given, it reorders data_df columns accordingly.

    :param data_df: a pandas dataframe with rate constants info
    :param column_order: column order in integers (optional)
    :param convert_to_ratios: whether or not the rate constants are converted into ratios.
    :return data_df: pandas dataframe with column names defined.
    """

    if convert_to_ratios:
        column_labels = [''.join(['$K_', str(i), '$'])  for i in range(1, len(data_df.columns) + 1)]
    else:
        column_labels = [''.join(['$\over', 'left' if i % 2 == 0 else 'right', 'arrow{k_{',  str(i // 2), '}}$']) for i in
                         range(2, len(data_df.columns)+2)]
    data_df.columns = column_labels

    if column_order:
        columns_reordered = [column_labels[entry] for entry in column_order]
        data_df = data_df[columns_reordered]

    #data_df.columns = column_labels[:2]

    return data_df


def import_rateconstants(file_in, ssd_threshold, Keq=None, column_order=None, convert_to_ratios=False, filter=True,
                         limit=500):
    """
    Import rate constant data from a tab separated file as a pandas dataframe where each row corresponds to a set of parameters and respective
    sum of squared deviations (ssd). The first column is the ssd and the remaining columns corresponds to rate constants,
    for a reaction with 3 steps, the expected order for the columns is the following: {ssd, k_1r, k_1f, k_2r, k_2f, k_3r, k_3f}.

    The rows with ssd < ssd_threshold are discarded, and the columns are reordered according to column_order, if specified.
    The data may also be converted to rate constants ratios of the form k_if/k_ir.
    The column will be named.


    :param file_in: path to tab separated file containing rate constant data
    :param ssd_threshold: ssd value below which rate constant sets will be discarded
    :param Keq:
    :param column_order: column order to reorder dataframe columns (in column indices)
    :param convert_to_ratios: whether or not the rate constants should be converted to rate constants ratios
    :param filter: boolean value: whether or not the models should be filtered by ssd
    :param limit: integer specifying how many models should be considered
    :return data_df: processed dataframe
    :return fitness: the list of of fitness values
    """
    data_df = pd.read_csv(file_in, header=None, sep='\t')
    fitness = data_df[0].values

    if filter:
        data_df = data_df[data_df[0] < ssd_threshold]
        fitness = data_df[0].values
    del data_df[0]

    # remove last columns for TalB, so that the Ki ratio is ignored
    if file_in.find('TAL') != -1:
        data_df = data_df[data_df.columns[:-2]]

    if convert_to_ratios:
        data_df = convert_df_to_ratios(data_df)


        n_cols = len(data_df.columns)
        n_keqs = len(column_order)
        for i in range(n_cols-n_keqs):
            del data_df[n_cols - i]

        if Keq:
            for i in data_df.index:
                if not (Keq * 0.9 < np.prod(data_df.ix[i, :]) < Keq * 1.1):
                    print i, len(fitness), np.prod(data_df.ix[i, :]), fitness[0], fitness[i]
                    print '----'
                    data_df = data_df.drop([i])

    if len(data_df.index) > limit:
        if limit < 100:
            ind_list = range(len(data_df.index))
            random.shuffle(ind_list)
            ind_list = ind_list[:limit]
            data_df = data_df.iloc[ind_list, :]
        else:
            data_df = data_df.iloc[:limit, :]
            fitness = fitness[:limit]

    data_df = pd.DataFrame(np.log10(data_df))

    data_df = set_df_column_labels(data_df, column_order=column_order, convert_to_ratios=convert_to_ratios)
    #print len(data_df.index)
    return data_df, fitness


def convert_df_to_ratios(data_df):
    """
    Converts rate constants data in pandas dataframe to ratios of rate constants, e.g. K_1 = k_1f/k_1r.

    :param data_df: pandas dataframe with rate constants data
    :return ratio_df: pandas dataframe with rate constants ratios.
    """

    ratio_df = pd.DataFrame()
    for i in range(1, len(data_df.columns)/2 + 1):
        ratio_df[i] = data_df[i*2]/data_df[i*2-1]

    return ratio_df


def import_data_timecourses(file_in, Keq=None, sub_list=None, prod_list=None, n_mets=5):
    """
    Given the path to the file with time-course data, imports data

    :param file_in: path to the file with time-course data
    :return: data frame with time-course data
    """

    data_df = pd.read_csv(file_in, header=None, sep='\t')

    data_df.columns = data_df.ix[0, :]
    data_df = data_df.drop(data_df.index[[0]])

    data_df = data_df.convert_objects(convert_numeric=True)

    #if Keq:
    #print 'Keq', Keq

    #last_t = len(data_df.index)

    #subs_product = np.prod([data_df.ix[last_t-1, sub].values for sub in sub_list], axis=0)
    #prods_product = np.prod([data_df.ix[last_t-1, prod].values for prod in prod_list], axis=0)
    #Keq_list = prods_product / subs_product

    #print Keq_list

    #models_to_remove = np.where((Keq_list > (Keq * 1.01)) | (Keq_list < (Keq * 0.99)))

    #print models_to_remove

    #assert len(models_to_remove[0]) == 0

    return data_df


# not in use
def analyze_running_times():
    file_in_nodKd = '/home/mrama/Dropbox/Kinetics/Enzymes_new/TALA2_dKd_statistics/results/raw/TALA2_nodKd_running_times.dat'
    file_in_dKd = '/home/mrama/Dropbox/Kinetics/Enzymes_new/TALA2_dKd_statistics/results/raw/TALA2_dKd_running_times.dat'
    file_out = '/home/mrama/Dropbox/Kinetics/Enzymes_new/TALA2_dKd_statistics/results/raw/TALA2_running_times_sumup.dat'

    data_df_nodKd = pd.read_csv(file_in_nodKd, header=None)
    data_df_dKd = pd.read_csv(file_in_dKd, header=None)

    n_fits = 4.
    with open(file_out, 'w+') as f_out:

        f_out.write('nodKd:\n')
        mean = (data_df_nodKd[0].mean() / n_fits) / 3600.
        std = (data_df_nodKd[0].std(ddof=1) / np.sqrt(n_fits)) / 3600.
        conf_int = stats.norm.interval(0.95, loc=mean, scale=std)

        f_out.write(''.join(['average running time per fit(h):\t%.2f' % mean, '\n']))
        f_out.write(''.join(['std running time per fit (percentage of mean):\t%.2f' % ((std / mean) * 100), '\n']))
        f_out.write(''.join(['95% confidence interval running time per fit (h):\t', str(conf_int), '\n']))


        f_out.write('dKd:\n')
        mean = (data_df_dKd[0].mean() / n_fits) / 3600.
        std = (data_df_dKd[0].std(ddof=1) / np.sqrt(n_fits)) / 3600.
        conf_int = stats.norm.interval(0.95, loc=mean, scale=std)

        f_out.write(''.join(['average running time per fit(h):\t%.2f' % mean, '\n']))
        f_out.write(''.join(['std running time per fit (percentage of mean):\t%.2f' % ((std / mean) * 100), '\n']))
        f_out.write(''.join(['95% confidence interval running time per fit (h):\t', str(conf_int), '\n']))

