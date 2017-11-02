import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from src.kinetics_integration.plots_definitions import ssd_heatmap_defs


def _import_data(file_in_list):

    data_df_median = pd.DataFrame()
    for file_in in file_in_list:
        data_df = pd.read_csv(file_in, sep=',', index_col=0, header=None)
        '''
        from astropy.stats import median_absolute_deviation
        print(data_df)
        print(data_df.transpose())

        print(data_df.transpose().median(axis=0))
        print(data_df.transpose().apply(median_absolute_deviation, axis=0).values / data_df.transpose().median(axis=0).values )

        print(data_df.transpose().apply(median_absolute_deviation, axis=0))
        print(np.max(data_df.transpose().apply(median_absolute_deviation, axis=0).values / data_df.transpose().median(axis=0).values )*100)
        exit()
        '''
        df_median = pd.DataFrame(data_df.transpose().median(axis=0))
        data_df_median = data_df_median.append(df_median.transpose())

    return data_df_median


def plot_ssd_heatmap(file_in_list, plot_file, enzyme_list):
    """
    Plot heatmap for sum of squared deviations.

    :param file_in_list: list of files (one per enzyme) that contain the ssd for each dKd value.
    :param plot_file: file path and name for the plot.
    :param enzyme_list: list of enzymes for whose ssd will be plotted.
    :return: None
    """

    ssd_heatmap_defs()

    data_df_median = _import_data(file_in_list)
    data_df_median.index = enzyme_list

    fig, ax = plt.subplots(1, 1, figsize=(8.5, 4))
    print [type(data_df_median.columns.values[i]) for i in range(len(data_df_median.columns.values))]

    print data_df_median.columns.values[1:]
    xtick_labels = list(np.log10([float(item) for item in data_df_median.columns.values[1:]]))
    xtick_labels = [str('$ 10^{' + str(item)[:-2] +'}$') for item in xtick_labels]
    xtick_labels.insert(0, 'ref')
    print xtick_labels

    cbar_ax = fig.add_axes([0.3, 0.95, 0.5, 0.04])

    sns.heatmap(data_df_median.apply(np.log10), ax=ax, xticklabels=xtick_labels, cbar_ax=cbar_ax,
                cmap='RdBu_r', center=0, vmin=-17, vmax=17, cbar_kws={'orientation': 'horizontal'})

    ax.set_xlabel('$\Delta K_b$ value')
    plt.sca(ax)

    plt.xticks(rotation=90)
    plt.tick_params(axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom='off',      # ticks along the bottom edge are off
                    top='off')       # labels along the bottom edge are off

    plt.yticks(rotation=0)
    plt.tick_params(axis='y',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    right='off',      # ticks along the bottom edge are off
                    left='off')       # labels along the bottom edge are off

    cbar_ax.tick_params(labelsize=20, length=0)

    fig.subplots_adjust(top=0.83, bottom=0.33, left=0.12, right=0.99, wspace=0.4)
    plt.savefig(''.join([plot_file, '.pdf']), dpi=300)
    plt.close()
    print(''.join([plot_file, '.pdf']))


def param_scan_heatmap(main_dir, enzyme_name_list):
    """
    Import data and call function to plot the ssd heatmap for all enzymes.

    :param main_dir: path to folder where you have all folders related to enzyme models
    :param enzyme_name_list: list with the names of enzymes for which you're doing this analysis.
    :return: None
    """

    file_in_list = []
    for enzyme_name in enzyme_name_list:
        file_in = ''.join([main_dir, enzyme_name, '/', enzyme_name, '_param_scan/output/treated_data/param_scan_ssd.csv'])
        file_in_list.append(file_in)

    plot_file = ''.join([main_dir, 'plots/param_scan_heatmap_ssd'])

    enzyme_list = ['ENO', 'GAPD', 'TALB']
    plot_ssd_heatmap(file_in_list, plot_file, enzyme_list)


def main():

    main_dir = '/home/mrama/Dropbox/PhD_stuff/Projects/MD/eMASS-MD/enzyme_models/'
    enzyme_name_list = ['ENO', 'GAPD', 'TALA2']
    param_scan_heatmap(main_dir, enzyme_name_list)


if __name__ == '__main__':
    main()

