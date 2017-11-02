import os

from src.kinetics_integration.analyze_var_IQR import calculate_elementary_keqs_iqr
from src.kinetics_integration.analyze_entropy import calculate_enzyme_entropy
from src.kinetics_integration.analyze_model_fitness import plot_fitness_dist, plot_n_good_models_dist
from src.kinetics_integration.analyze_time_courses import time_course_analysis
from src.kinetics_integration.plot_clustermaps import cluster_map_param_inf


def run_fitness_analysis(base_folder, file_in_base, model_type_list, ssd_threshold, Keq, column_order,
                         column_labels, x_label, scale_data, convert_to_ratios, n_model_sets, param_inf_type):
    color_list = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']

    if not os.path.exists(''.join([base_folder, 'plots'])):
        os.makedirs(''.join([base_folder, 'plots']))

    Keq_local = None
    #print 'in', Keq_local
    file_out_base = ''.join([base_folder, 'plots/ENO_fitness_dist_all_', param_inf_type])
    plot_fitness_dist(file_in_base, file_out_base, model_type_list, n_model_sets, ssd_threshold, Keq_local,
                      column_order, column_labels, x_label, convert_to_ratios, color_list, filter=False)

    Keq_local = Keq
    #print 'in', Keq_local
    file_out_base = ''.join([base_folder, 'plots/ENO_n_valid_models_Keq_check_', param_inf_type])
    plot_n_good_models_dist(file_in_base, file_out_base, model_type_list, n_model_sets, ssd_threshold, Keq_local,
                            column_order,
                            column_labels, x_label, scale_data, convert_to_ratios, color_list, filter=True)

    Keq_local = None
    #print 'in', Keq_local
    file_out_base = ''.join([base_folder, 'plots/ENO_n_valid_models_NO_Keq_check_', param_inf_type])
    plot_n_good_models_dist(file_in_base, file_out_base, model_type_list, n_model_sets, ssd_threshold, Keq_local,
                            column_order,
                            column_labels, x_label, scale_data, convert_to_ratios, color_list, filter=True)


def analyze_eno(base_folder, enzyme, n_model_ensembles, Keq, e_total, clustermaps, fitness_analysis, keq_range_analysis,
                entropy_analysis, time_courses, param_inf_type):

    ssd_threshold = 1
    scale_data = False

    column_order = [0, 2, 1]

    n_variables = len(column_order)
    convert_to_ratios = True

    model_type_list = ['all', 'dKd', 'Keq', 'Km', 'kcat']
    column_labels = ['None', '$\Delta K_d$', '$K_{eq}$', '$K_m^{2pg}$', '$k_{cat}$']
    x_label = 'Data point removed'

    file_in_base = ''.join([base_folder, 'treated_data/rateconst_ENO_'])

    if clustermaps:
        print 'plotting clustermaps'

        if not os.path.exists(''.join([base_folder, 'clustermaps'])):
            os.makedirs(''.join([base_folder, 'clustermaps']))

        vmin = -6
        vmax = 9

        model_types = ['all', 'dKd']
        #n_model_ensembles = 100

        Keq_local = None
        fig_size = (4, 6)
        file_out_base = ''.join([base_folder, 'clustermaps/ENO_rateconst_', param_inf_type])
        column_order_rateconst = [0, 1, 4, 5, 2, 3]
        cluster_map_param_inf(file_in_base, file_out_base, ssd_threshold, False, column_order_rateconst, vmin, vmax,
                              model_types, n_model_ensembles, fig_size, Keq=Keq_local)

    if fitness_analysis:
        print 'doing fitness analysis'

        run_fitness_analysis(base_folder, file_in_base, model_type_list, ssd_threshold, Keq, column_order,
                             column_labels, x_label, scale_data, convert_to_ratios, n_model_ensembles, param_inf_type)

    if keq_range_analysis:
        print 'analyzing variables IQR'

        if param_inf_type == 'Keq':
            zorder_list = [2, 0, 1, 3]
        else:
            zorder_list = [4, 3, 0, 1, 2]


        color_list = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']

        legend_details = [(0.5, 1.15), 0.875]
        legend_n_cols = 5

        Keq_local = None
        n_model_ensembles = 100

        file_out_base = ''.join([base_folder, 'plots/ENO_elementary_rateconst_range_', param_inf_type])
        column_order_rateconst = [0, 1, 4, 5, 2, 3]

        x_labels = ['$\overleftarrow{k_{1}}$', '$\overrightarrow{k_{1}}$', '$\overleftarrow{k_{2}}$',
                    '$\overrightarrow{k_{2}}$', '$\overleftarrow{k_{3}}$', '$\overrightarrow{k_{3}}$']

        calculate_elementary_keqs_iqr(file_in_base, file_out_base, model_type_list, n_model_ensembles, ssd_threshold,
                                      Keq_local, column_order_rateconst, False, column_labels, x_labels, color_list,
                                      zorder_list)

    if entropy_analysis:
        print 'analyzing entropy'

        if not os.path.exists(''.join([base_folder, 'entropy'])):
            os.makedirs(''.join([base_folder, 'entropy']))

        color_list = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']

        plot_entropy_only = False
        bin_width_list = [3]
        fixed_size = False
        y_lims = [0.1, 0.5]
        Keq_local = None

        # start try
        #convert_to_ratios = False
        #column_order = [0, 2, 1]
        #column_order_rateconst = [0, 1, 4, 5, 2, 3]

        # end try
        """
        file_out_base = ''.join([base_folder, 'entropy/ENO_ratios_'])
        limit = 100
        calculate_enzyme_entropy(file_in_base, file_out_base, n_model_ensembles, model_type_list, column_order,
                                 n_variables, ssd_threshold, Keq_local, convert_to_ratios, limit, bin_width_list,
                                 column_labels, y_lims, x_label, color_list, plot_entropy_only=plot_entropy_only,
                                 fixed_size=fixed_size, n_samples_per_bin=None)

        limit = 45

        calculate_enzyme_entropy(file_in_base, file_out_base, n_model_ensembles, model_type_list, column_order,
                                 n_variables, ssd_threshold, Keq_local, convert_to_ratios, limit, bin_width_list,
                                 column_labels, y_lims, x_label, color_list, plot_entropy_only=plot_entropy_only,
                                 fixed_size=fixed_size, n_samples_per_bin=None)
        """

        file_out_base = ''.join([base_folder, 'entropy/ENO_rateconst_'])
        column_order_rateconst = [0, 1, 4, 5, 2, 3]
        n_variables_rateconst = len(column_order_rateconst)
        limit = 100
        calculate_enzyme_entropy(file_in_base, file_out_base, n_model_ensembles, model_type_list,
                                 column_order_rateconst, n_variables_rateconst, ssd_threshold, Keq_local, False, limit,
                                 bin_width_list, column_labels, y_lims, x_label, color_list,
                                 plot_entropy_only=plot_entropy_only, fixed_size=fixed_size, n_samples_per_bin=None)

        limit = 45

        calculate_enzyme_entropy(file_in_base, file_out_base, n_model_ensembles, model_type_list,
                                 column_order_rateconst, n_variables_rateconst, ssd_threshold, Keq_local, False, limit,
                                 bin_width_list, column_labels, y_lims, x_label, color_list,
                                 plot_entropy_only=plot_entropy_only, fixed_size=fixed_size, n_samples_per_bin=None)

    if time_courses:
        print 'analyzing time-courses'
        plot_only = False

        subs_list = ['2pg']
        prod_list = ['pep']
        perturbation_list = ['noPerturb']

        filter_models = False
        x_lims = [0, 10**3]

        time_point_list = [0, 100]
        plot_type_list = ['mets', 'enz', 'flux']

        for perturb in perturbation_list:
            print perturb
            for plot_type in plot_type_list:
                print plot_type
                mets_order = None
                legend_label = None

                if plot_type == 'mets':

                    file_type = 'conc_mets'
                    #y_lims = [0, 0.0025]  # 0.00022
                    species_list = ['2pg', 'pep']
                    mets_order = species_list
                    y_label = 'Concentration (mol/L)'
                    y_lims = [(0, 0.0001), (0.00015, 0.00025)]
                    plot_spacings = [0.89, 0.15, 0.15, 0.98]
                    fig_size = (6, 4)

                elif plot_type == 'enz':
                    file_type = 'conc_enz'
                    #y_lims = [0, 7 * 10 ** -5]
                    species_list = ['ENO', 'ENO&2pg', 'ENO&pep']
                    legend_label = ['E', 'E-2pg', 'E-pep']
                    y_label = 'Concentration (mol/L)'
                    y_lims = [[0, 0.00004], [0, 0.000008], [0, 0.00008], [0., 0.00008], [0., 0.00003], [0., 0.00003]]
                    plot_spacings = [0.89, 0.2, 0.15, 0.98]
                    fig_size = (6, 4)

                elif plot_type == 'flux':
                    file_type = 'flux'
                    #y_lims = [-0.001, 0.001]
                    species_list = ['v_eno1', 'v_eno2', 'v_eno3']
                    legend_label = ['eno1', 'eno2', 'eno3']
                    y_label = 'Flux (mol/L/s)'
                    y_lims = None
                    plot_spacings = [0.89, 0.24, 0.15, 0.98]
                    fig_size = (6, 4)

                else:
                    print 'plot_type not specified correctly, ', plot_type
                    return

                file_in_base = ''.join([base_folder, 'model_simulations/data/', perturb, '/sim_res_', file_type])
                file_out_base = ''.join([base_folder, 'model_simulations/plots/', perturb, '/', plot_type, '/sim_res_', file_type])

                if not os.path.exists(''.join([base_folder, 'model_simulations/plots/', perturb, '/', plot_type])):
                    os.makedirs(''.join([base_folder, 'model_simulations/plots/', perturb, '/', plot_type]))

                time_course_analysis(file_in_base, file_out_base, species_list, Keq, subs_list, prod_list, mets_order,
                                     legend_label, n_model_ensembles, plot_type, time_point_list, plot_spacings,
                                     fig_size, y_lims=y_lims, y_lims_timecourses=y_lims, plot_only=plot_only)
if __name__ == '__main__':
    enzyme = 'ENO'
    e_total = 1.93*10**-5
    Keq = 5.19
    #Keq = ''
    param_inf_type = ''

    clustermaps = True
    fitness_analysis = True
    keq_range_analysis = True
    entropy_analysis = True
    time_courses = False
    n_model_ensembles = 100

    base_folder = '/home/mrama/Dropbox/PhD_stuff/Projects/MD/eMASS-MD/enzyme_models/ENO/ENO_param_inf/output/'

    analyze_eno(base_folder, enzyme, n_model_ensembles, Keq, e_total, clustermaps, fitness_analysis, keq_range_analysis,
                entropy_analysis, time_courses, param_inf_type)
