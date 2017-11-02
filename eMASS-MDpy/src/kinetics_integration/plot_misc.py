"""
Not in use
"""

from collections import OrderedDict
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import ticker as tck
from statsmodels.robust.scale import mad
from math import log10


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold')
plt.rcParams['axes.linewidth'] = 1.5

plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['xtick.minor.width'] = 1.5
plt.rcParams['ytick.minor.width'] = 1.0

plt.rcParams['xtick.major.pad'] = '12'
plt.rcParams['ytick.major.pad'] = '8'


def plot_dKd_prediction(file_in, plot_file, true_dKd, Keq, plot_type='hist'):

    axes_label_size = 14

    data_df = pd.read_csv(file_in, sep='\t', header=None)
    print data_df
    #print 10**np.log10(data_df[0]).mean(), data_df[0].mean()
    #print 10**np.log10(data_df[0]).median(), data_df[0].median()

    #print data_df[(data_df[0] > 0.2) & (data_df[0] < 0.3)].count() / 10000
    #print data_df[(data_df[0] > 0.) & (data_df[0] < 1)].count() / 10000

    if plot_type == 'hist':
        ax = np.log10(data_df[0]).hist(bins=100)
        plt.plot([true_dKd, true_dKd], [0, 400], 'k-', color='red', lw=4)
        plt.plot([Keq, Keq], [0, 400], 'k-', color='green', lw=4)
        plt.xlabel('$\Delta K_d$ order of magnitude ($10^x$)', size=axes_label_size, labelpad=5)

    elif plot_type == 'boxplot':
        ax = sns.boxplot(data=np.log10(data_df[0]), linewidth=1.7, saturation=1, flierprops=dict(linestyle='none'))
        plt.plot([-1, 1], [true_dKd, true_dKd], 'k-', color='red', lw=4)
        plt.plot([-1, 1], [Keq, Keq], 'k-', color='green', lw=4)
        plt.ylabel('$\Delta K_d$ order of magnitude ($10^x$)', size=axes_label_size, labelpad=5)

    plt.tight_layout()
    plt.savefig(''.join([plot_file, '_', plot_type, '.pdf']), dpi=300)
    plt.savefig(''.join([plot_file, '_', plot_type, '.png']), dpi=300)
    plt.close()


def plot_eigenvalues_dist(file_in, plot_file):

    axes_label_size = 14

    data_df = pd.read_csv(file_in, sep='\t', header=None)
    ax = data_df[0].hist(bins=100)
    ax.set_xlim(-2 * 10**9, 0.2 * 10**9)

    print data_df[(data_df[0] > 0)].count() / 10000
    print data_df[(data_df[0] > 0)].max()

    plt.xlabel('eigenvalues', size=axes_label_size, labelpad=5)
    plt.tight_layout()
    plt.savefig(''.join([plot_file, '.pdf']), dpi=300)
    plt.savefig(''.join([plot_file, '.png']), dpi=300)
    plt.close()


def ENO_plots():
    enzyme = 'ENO'
    true_dKd = np.log10(0.25)
    Keq = np.log10(5.19)
    nodKd_fit = 3.3 * 10 ** -22

    file_in = '/home/mrama/Dropbox/Kinetics/Enzymes_new/ENO_dKd_statistics/results/treated_data/dKdPredicted_ENO_nodKd_all.csv'
    plot_file = '/home/mrama/Desktop/meeting_july21/dKdPredicted_ENO'
    plot_dKd_prediction(file_in, plot_file, true_dKd, Keq, plot_type='hist')
    plot_dKd_prediction(file_in, plot_file, true_dKd, Keq, plot_type='boxplot')

    #file_in = '/home/mrama/Dropbox/Kinetics/Enzymes_new/ENO_dKd_statistics/results/treated_data/eigenvalues_dKd_singleSet.csv'
    #plot_file = '/home/mrama/Desktop/meeting_july8/eigenvalues_ENO_dKd'
    #plot_eigenvalues_dist(file_in, plot_file)

    #file_in = '/home/mrama/Dropbox/Kinetics/Enzymes_new/ENO_dKd_statistics/results/treated_data/eigenvalues_NO_dKd_singleSet.csv'
    #plot_file = '/home/mrama/Desktop/meeting_july8/eigenvalues_ENO_NO_dKd'
    #plot_eigenvalues_dist(file_in, plot_file)

def G6PD_plots():
    enzyme = 'G6PDH2r'
    true_dKd = 23
    Keq = np.log10(18.3)
    nodKd_fit = 3.3 * 10 ** -22

    file_in = '/home/mrama/Dropbox/Kinetics/Enzymes_new/G6PDH2r_dKd_statistics/results/treated_data/dKdPredicted_G6PDH2r_nodKd_all.csv'
    plot_file = '/home/mrama/Desktop/meeting_july21/dKdPredicted_G6PDH2r'
    plot_dKd_prediction(file_in, plot_file, true_dKd, Keq, plot_type='hist')
    plot_dKd_prediction(file_in, plot_file, true_dKd, Keq, plot_type='boxplot')


def GAPD_plots():
    enzyme = 'GAPD'
    true_dKd = np.log10(7.53 * 10**4)
    Keq = np.log10(0.408)
    #nodKd_fit = 3.3 * 10 ** -22

    file_in = '/home/mrama/Dropbox/Kinetics/Enzymes_new/GAPD_dKd_statistics/results/treated_data/dKdPredicted_GAPD_nodKd_all.csv'
    plot_file = '/home/mrama/Desktop/meeting_july21/dKdPredicted_GAPD'
    plot_dKd_prediction(file_in, plot_file, true_dKd, Keq, plot_type='hist')
    plot_dKd_prediction(file_in, plot_file, true_dKd, Keq, plot_type='boxplot')


def TALB_plots():
    enzyme = 'TALA2'
    true_dKd = np.log10(125.8)
    Keq = np.log10(0.84)
    #nodKd_fit = 3.3 * 10 ** -22

    file_in = '/home/mrama/Dropbox/Kinetics/Enzymes_new/TALA2_dKd_statistics/results/treated_data/dKdPredicted_TALA2_nodKd_all.csv'
    plot_file = '/home/mrama/Desktop/meeting_july21/dKdPredicted_TALA2'
    plot_dKd_prediction(file_in, plot_file, true_dKd, Keq, plot_type='hist')
    plot_dKd_prediction(file_in, plot_file, true_dKd, Keq, plot_type='boxplot')

if __name__ == '__main__':
    #ENO_plots()
    #G6PD_plots()
    #GAPD_plots()
    TALB_plots()
