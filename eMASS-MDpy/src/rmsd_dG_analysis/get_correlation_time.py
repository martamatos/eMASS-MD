import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.tools.plotting import autocorrelation_plot

def _get_correlation(mean_list):

    mean_list_mean = np.mean(mean_list)
    corr = np.sum([((np.mean(mean_list[0] * dg_val) - mean_list_mean ** 2) / (np.mean(np.array(mean_list) ** 2) - mean_list_mean ** 2)) for dg_val in mean_list])

    return(corr)


def _get_block_var(dg_data_df, block_size):
    mean_list = []
    for i in range(0, len(dg_data_df.index), block_size):
        # print(dg_data_df['deltaG'][i:i+block_size])
        # print(dg_data_df['deltaG'][i:i+block_size].mean())
        mean_list.append(dg_data_df['deltaG'][i:i + block_size].mean())
        # if block_size == 70:
        # print('*************')
        # print(dg_data_df['deltaG'][i:i+block_size])
        # print(dg_data_df['deltaG'][i:i+block_size].mean())
        # print(dg_data_df['deltaG'][i:i+block_size].var())
        # print('*************')

    corr = _get_correlation(mean_list)
    n_samples = len(mean_list)
    dg_mean = np.mean(mean_list)
    dg_values = mean_list
    var = (1 /(n_samples * (n_samples - 1.0))) * np.sum([(dg_val - dg_mean)**2 for dg_val in dg_values])

    return(var, corr)


def get_statistical_inefficiency(dg_data_df, block_size_list):
    # dg_var = dg_data_df['deltaG'].var()
    dg_values = dg_data_df['deltaG'].values
    dg_mean = dg_data_df['deltaG'].mean()
    n_samples = 70
    dg_var = (1 /(n_samples * (n_samples - 1.0))) * np.sum([(dg_val - dg_mean)**2 for dg_val in dg_values])

    print(dg_data_df['deltaG'].mean())
    print(dg_data_df['deltaG'].std())
    print(dg_var)
    print('-------')

    stat_inefficiency_list = []
    corr_list = []
    for block_size in block_size_list:
        (block_var, corr) = _get_block_var(dg_data_df, block_size)
        stat_inefficiency_list.append((block_var / dg_var))
        corr_list.append(corr)

    print(stat_inefficiency_list)
    print(corr_list)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].scatter(x=block_size_list, y=stat_inefficiency_list)
    ax[0].set_title('stat inef')
    ax[0].set_xlim(-1, 36)
    ax[0].set_ylim(-10, 10)

    ax[1].scatter(x=block_size_list, y=corr_list)
    ax[1].set_title('corr')
    ax[1].set_xlim(-1, 71)
    ax[1].set_ylim(-10, 36)
    plt.show()
    plt.close()

    return 0

dg_data_df = pd.read_csv(
    '/home/mrama/Desktop/MD/eMASS-MD_complete_data/MD_data/GAPDH/4_MMPBSA2016/1_NAD/1_G3P/GAPDH_WT_NAD_G3P_docked_cl000.dat',
    sep='\t', index_col='frame')
print(dg_data_df.ix[:70])

block_size_list = range(1, 35)

#get_statistical_inefficiency(dg_data_df.ix[10:70], block_size_list)

import random
data = [random.random() for i in range(70)]

from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.tsa.stattools import acf

print(acf(data))

fig, ax = plt.subplots(2, sharex=True, sharey=True)
plot_acf(dg_data_df.ix[8:70:1], ax=ax[0])


plot_acf(dg_data_df.ix[8:70:2], ax=ax[1])
plt.show()

exit()

def acorr(x, ax=None):
    if ax is None:
        ax = plt.gca()

    #x = x - np.mean(x)

    autocorr = np.correlate(x, x, mode='same')
    #autocorr = autocorr[x.size:]
    #autocorr /= autocorr.max()
    print(autocorr)

    return ax.scatter(x=range(len(autocorr)), y=autocorr)


#acorr(data)
autocorrelation_plot(data)
acorr(data)
plt.show()
exit()
"""
fig, ax = plt.subplots(9, 1, figsize=(20, 17), sharex=True, sharey=True)
autocorrelation_plot(dg_data_df.ix[0:70:1], ax=ax[0])
autocorrelation_plot(dg_data_df.ix[8:70:1], ax=ax[1])
autocorrelation_plot(dg_data_df.ix[8:70:2], ax=ax[2])
autocorrelation_plot(dg_data_df.ix[8:70:3], ax=ax[3])
autocorrelation_plot(dg_data_df.ix[8:70:4], ax=ax[4])
autocorrelation_plot(dg_data_df.ix[8:70:5], ax=ax[5])
autocorrelation_plot(dg_data_df.ix[8:70:6], ax=ax[6])
autocorrelation_plot(dg_data_df.ix[8:70:7], ax=ax[7])
autocorrelation_plot(dg_data_df.ix[8:70:8], ax=ax[8])
plt.xticks(range(70), rotation=90)
plt.savefig('/home/mrama/PhD_stuff/Projects/MD/autocorrelation_diff_space_start8.pdf')
plt.tight_layout()
plt.close()"""

data = [item[0] for item in dg_data_df.ix[8:70:2].values]

from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.tsa.stattools import acf

print(acf(data))

fig, ax = plt.subplots(3, sharex=True, sharey=True)
autocorrelation_plot(dg_data_df.ix[8:70:2], ax=ax[0])


res = plot_acf(data, ax=ax[1])
print(res)




def acorr(x, ax=None):
    if ax is None:
        ax = plt.gca()

    x = x - np.mean(x)

    autocorr = np.correlate(x, x, mode='valid')
    #autocorr = autocorr[x.size:]
    #autocorr /= autocorr.max()
    print(autocorr)

    return ax.scatter(x=range(len(autocorr)), y=autocorr)


acorr(data, ax[2])
ax[0].set_xlim(0,32)
plt.show()

plt.close()



