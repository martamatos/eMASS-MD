import pandas as pd
import seaborn as sns

COLOR_LIST = sns.color_palette("Set1")


def get_keq_bins_fixed_width(data_df, bin_width):
    """
    Given a pandas dataframe and a bin width, it bins the dataframe columns using bins with the width specified.

    :param data_df: pandas dataframe
    :param bin_width: an integer specifying the bin's width in orders of magnitude
    :return pattern_count_dic: a dictionary where the index are the different patterns and the first column is the
                                frequency
    :return  pattern_dic_by_model: a dictionary where the columns are the patterns and the values are the models that
                                   have that specific pattern
    """

    n_keqs = len(data_df.columns)
    pattern_count_dic = {}
    pattern_dic_by_model = {}
    max_keq = data_df.max().max()
    min_keq = abs(data_df.min().min())
    assert 10 ** min_keq >= 10. ** -15 and 10 ** max_keq <= 10. ** 15

    n_pos_bins = int(max_keq / bin_width) if max_keq % bin_width == 0 else int(max_keq / bin_width) + 1
    n_neg_bins = int(min_keq / bin_width) if min_keq % bin_width == 0 else int(min_keq / bin_width) + 1

    pos_bins_limits = [[i * bin_width, (i + 1) * bin_width] for i in range(n_pos_bins)]

    neg_bins_limits = [[-i * bin_width, -(i + 1) * bin_width] for i in range(n_neg_bins)]

    for i in data_df.index.values:
        pattern_temp = []
        for j in range(n_keqs):
            keq_val = data_df.ix[i, j]
            if keq_val < 0:
                for bin_i in range(n_neg_bins):
                    if neg_bins_limits[bin_i][0] >= keq_val > neg_bins_limits[bin_i][1]:
                        pattern_temp.append(-bin_i - 1)
            elif keq_val >= 0:
                for bin_i in range(n_pos_bins):
                    if pos_bins_limits[bin_i][0] <= keq_val < pos_bins_limits[bin_i][1]:
                        pattern_temp.append(bin_i + 1)
            else:
                'oooops------------'
                'oooops------------'
                'oooops------------'
                'oooops------------'
                'oooops------------'

        pattern_temp = tuple(pattern_temp)
        if pattern_temp in pattern_count_dic:
            pattern_count_dic[pattern_temp] += 1
            pattern_dic_by_model[pattern_temp].append(i)
        else:
            pattern_count_dic[pattern_temp] = 1
            pattern_dic_by_model[pattern_temp] = [i]

    return pattern_count_dic, pattern_dic_by_model


def _find_sign_change_entry(sorted_column, left, right):
    entry_gap = abs(right - left)

    while entry_gap > 1:

        if (sorted_column[left] < 0 and sorted_column[right] > 0) and entry_gap > 1:
            right = left + (right - left) / 2

        elif (sorted_column[left] > 0 and sorted_column[right] > 0) and entry_gap > 1:
            right_temp = left
            left = left - (right - left)
            right = right_temp

        elif (sorted_column[left] < 0 and sorted_column[right] < 0) and entry_gap > 1:
            left_temp = right
            right = right + (right - left)
            right = len(sorted_column) - 1 if right >= len(sorted_column) else right
            left = left_temp

        entry_gap = abs(right - left)

    while sorted_column[left] < 0 and sorted_column[right] < 0:
        right += 1
    while sorted_column[left] > 0 and sorted_column[right] > 0:
        right -= 1

    return right


def _bin_fixed_size(bin_i, n_bins, n_samples, n_samples_per_bin, model_pattern_count_dic, sorted_column_index,
                    sorted_column):

    initial_bin_i = bin_i
    last_bin_ind = 0
    ind = 0
    samples_per_bin_count = None
    bin_borders = [sorted_column[0]]

    while bin_i < n_bins and ind < n_samples:
        interval_gap = 0
        samples_per_bin_count = 0

        while (samples_per_bin_count < n_samples_per_bin or interval_gap < 1) and ind < n_samples:
            model_pattern_count_dic[sorted_column_index[ind]].append(bin_i)
            interval_gap = sorted_column[ind] - sorted_column[last_bin_ind]
            ind += 1
            samples_per_bin_count += 1

        last_bin_ind = ind - 1
        bin_borders.append(sorted_column[last_bin_ind])
        bin_i += 1

    if samples_per_bin_count and samples_per_bin_count < n_samples_per_bin / 2 + 1 and bin_i > initial_bin_i + 1:
        for i in range(1, samples_per_bin_count + 1):
            ind = -i
            model_pattern_count_dic[sorted_column_index[ind]][-1] = bin_i - 2

        bin_i -= 1
        bin_borders[-2] = sorted_column[-1]
        bin_borders = bin_borders[:-1]

    return bin_i, bin_borders


def get_keq_bins_fixed_size(data_df, file_out, n_samples_per_bin=5):
    """
    Given a pandas dataframe and the number of samples that each bin should contain, it bins the dataframe's columns
    using bins with the number of samples specified.

    :param data_df: pandas dataframe
    :param file_out: file base name where binning results (also intermediate) will be stored.
    :param n_samples_per_bin: an integer how many samples should be in each bin
    :return pattern_count_dic: a dictionary where the index are the different patterns and the first column is the
                                frequency
    :return  pattern_dic_by_model: a dictionary where the columns are the patterns and the values are the models that
                                    have that specific pattern
    """

    n_keqs = len(data_df.columns)
    n_samples = len(data_df.index)
    remainder = n_samples % n_samples_per_bin
    n_bins = n_samples / n_samples_per_bin + 1 if remainder != 0 else n_samples / n_samples_per_bin

    bin_borders_dic = {i: [] for i in range(n_keqs)}
    model_pattern_count_dic = {sample_i: [] for sample_i in range(n_samples)}

    for i in range(n_keqs):
        sorted_column = data_df.ix[:, i].sort_values().values
        sorted_column_index = data_df.ix[:, i].sort_values().index

        if (sorted_column[0] > 0 and sorted_column[-1] < 0) or (sorted_column[0] < 0 and sorted_column[-1] > 0):
            left = 0
            right = n_samples / 2
            sign_change_entry = _find_sign_change_entry(sorted_column, left, right)
            assert sorted_column[sign_change_entry - 1] < 0 and sorted_column[sign_change_entry] > 0

            bin_i = 0
            bin_i, bin_borders_neg = _bin_fixed_size(bin_i, n_bins, sign_change_entry, n_samples_per_bin,
                                                     model_pattern_count_dic, sorted_column_index[:sign_change_entry],
                                                     sorted_column[:sign_change_entry])
            bin_i, bin_borders_pos = _bin_fixed_size(bin_i, n_bins, (n_samples - sign_change_entry), n_samples_per_bin,
                                                     model_pattern_count_dic, sorted_column_index[sign_change_entry:],
                                                     sorted_column[sign_change_entry:])

            bin_borders_dic[i] = bin_borders_neg + bin_borders_pos
        else:
            bin_i = 0
            bin_i, bin_borders_dic[i] = _bin_fixed_size(bin_i, n_bins, n_samples, n_samples_per_bin,
                                                        model_pattern_count_dic, sorted_column_index, sorted_column)

    assert len(model_pattern_count_dic.keys()) == n_samples
    bin_borders_df = pd.DataFrame().from_dict(bin_borders_dic, orient='index')
    bin_borders_df = bin_borders_df.transpose()
    bin_borders_df.to_csv(''.join([file_out, '_bin_borders.csv']))

    pattern_count_dic = {}
    pattern_dic_by_model = {}

    for i in range(n_samples):
        pattern = tuple(model_pattern_count_dic[i])
        if pattern in pattern_count_dic:
            pattern_dic_by_model[pattern].append(i)
            pattern_count_dic[pattern] += 1
        else:
            pattern_dic_by_model[pattern] = [i]
            pattern_count_dic[pattern] = 1

    return pattern_count_dic, pattern_dic_by_model


def bin_keqs(data_df, bin_width=None, fixed_size=False, n_samples_per_bin=5, file_out=None):
    """

    Given the bin_width and a dataframe with n sets of k elementary equilibrium constants  or rate constants (a model),
    it attributes each k_i constant to a bin, and each set of  constants
    is mapped to pattern (a set of k bins). In the end the n sets of k elementary equilibrium constants are all mapped
    to patterns, making it possible to calculate a probability distribution for each pattern.

    :param data_df: a pandas dataframe with a sets of elementary rate or equilibrium constants.
    :param bin_width: width of each bin in orders of magnitude
    :param fixed_size: boolean specifying whether the bins should have a fixed size, i.e. a fixed number of samples per bin.
    :param n_samples_per_bin: if fixed_size=True, how many samples should each bin have.
    :param file_out: file base name where binning results (also intermediate) will be stored if fixed_size=true.
    :return pattern_count_df: a dataframe with the number of models per pattern
    :return pattern_by_model_df: a dataframe mapping each set of elementary equilibrium constants to the respective
                                pattern
    """

    if fixed_size:
        pattern_count_dic, pattern_dic_by_model = get_keq_bins_fixed_size(data_df, file_out,
                                                                          n_samples_per_bin=n_samples_per_bin)
    elif bin_width:
        pattern_count_dic, pattern_dic_by_model = get_keq_bins_fixed_width(data_df, bin_width)

    else:
        print 'hummm... something went wrong on choosing the binning'
        return

    pattern_count_df = pd.DataFrame.from_dict(pattern_count_dic, orient='index')
    pattern_by_model_df = pd.DataFrame.from_dict(pattern_dic_by_model, orient='index')

    return pattern_count_df, pattern_by_model_df.transpose()
