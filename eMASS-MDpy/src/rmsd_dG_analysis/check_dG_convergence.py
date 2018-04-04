import pandas as pd
import random
import numpy as np
import matplotlib.pyplot as plt
from src.kinetics_integration.plots_definitions import get_elementary_keq_range_scatter_defs

COLOR_LIST = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']



def analyze_convergence(file_in, ligand_name, file_out, y_lims):

    get_elementary_keq_range_scatter_defs()

    data_df = pd.read_csv(file_in, sep='\t', index_col=0)
    data_df[data_df['mean'] > 0] = np.NaN

    dG_values = data_df.filter(regex=ligand_name, axis=0)['mean'].dropna().values
    n_dG_values = len(dG_values)
    ind_list = range(n_dG_values)

    fig, ax = plt.subplots()

    print(np.mean(dG_values))
    print('***********')
    n_samples = 1000
    #n_samples = n_dG_values
    sterr_last_time_point_list = []
    for i in range(1000):
        #random.shuffle(ind_list)
        ind_list = [random.randint(0, n_dG_values-1) for l in range(n_samples)]

        dG_values_random = dG_values[ind_list]
        #cumulative_mean = [np.mean(dG_values_random[:k]) for k in range(1, n_samples+1)]
        cumulative_sterr = [np.sqrt(np.var(dG_values_random[:k]))/float(np.sqrt(len(dG_values_random[:k]))) for k in range(1, n_samples+1)]
        sterr_last_time_point_list.append(cumulative_sterr[-1])
        ax.scatter(x=range(1, n_samples+1), y=cumulative_sterr, color=COLOR_LIST[1])

    print(min(sterr_last_time_point_list), np.mean(sterr_last_time_point_list), np.median(sterr_last_time_point_list), max(sterr_last_time_point_list))
    ax.set_xlim(0, n_samples+1)
    ax.set_ylim(y_lims)
    ax.set_xlabel('number of samples included in the mean value')
    ax.set_ylabel('mean dG kcal/mol')
    ax.grid(True)
    plt.tight_layout()

    if n_samples != n_dG_values:
        plt.savefig(''.join([file_out, '_', ligand_name, '_randint_1000_', str(n_samples), 'samples.png']))
    else:
        plt.savefig(''.join([file_out, '_', ligand_name, '.png']))
    plt.close()

    # len(cumulative_mean)
    #print cumulative_mean

    #_get_cumulative_mean(dG_values[random_ind])

    #print(data_df)
    return 0


def analyze_ENO(base_dir):

    file_in = ''.join([base_dir, 'ENO_AB/4_MMPBSA2016/ENO_ddG_sumup_1ns.csv'])
    file_out = ''.join([base_dir, 'ENO_AB/4_MMPBSA2016/ENO_convergence_stderr'])

    ligand_name = '_2PG_'
    y_lims = [0, 15]
    analyze_convergence(file_in, ligand_name, file_out, y_lims)

    ligand_name = '_PEP_'
    y_lims = [0, 15]
    analyze_convergence(file_in, ligand_name, file_out, y_lims)


def analyze_GAPD(base_dir):

    file_in = ''.join([base_dir, 'GAPDH/4_MMPBSA2016/GAPD_ddG_sumup_1ns.csv'])
    file_out = ''.join([base_dir, 'GAPDH/4_MMPBSA2016/GAPD_convergence_stderr'])

    ligand_name = '_G3P_'
    y_lims = [0, 15]
    analyze_convergence(file_in, ligand_name, file_out, y_lims)

    ligand_name = '_DPG_'
    y_lims = [0, 15]
    analyze_convergence(file_in, ligand_name, file_out, y_lims)


def analyze_TALB(base_dir):

    file_in = ''.join([base_dir, 'TalB/4_MMPBSA2016/TALB_ddG_sumup_1ns.csv'])
    file_out = ''.join([base_dir, 'TalB/4_MMPBSA2016/TALB_convergence_stderr'])

    ligand_name = '_S7P_linear'
    y_lims = [0, 15]
    analyze_convergence(file_in, ligand_name, file_out, y_lims)

    ligand_name = '_F6P_'
    y_lims = [0, 15]
    analyze_convergence(file_in, ligand_name, file_out, y_lims)


def analyze_TALB_old(base_dir):

    file_in = ''.join([base_dir, 'TalB/4_MMPBSA2016_old/TALB_ddG_sumup_20ns.csv'])
    file_out = ''.join([base_dir, 'TalB/4_MMPBSA2016_old/TALB_convergence_stderr'])

    ligand_name = '_S7P_linear'
    y_lims = [-40, 0]
    analyze_convergence(file_in, ligand_name, file_out, y_lims)

    ligand_name = '_F6P_'
    y_lims = [-40, 0]
    analyze_convergence(file_in, ligand_name, file_out, y_lims)




def main():

    base_dir = '/home/mrama/Desktop/MD/eMASS-MD_complete_data/MD_data/'
    analyze_TALB(base_dir)
    #analyze_TALB_old(base_dir)

if __name__ == '__main__':
    main()