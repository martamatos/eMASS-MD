from src.kinetics_integration.analyze_model_fitness import plot_fitness_dist, plot_n_good_models_dist
from src.kinetics_integration.analyze_var_IQR import calculate_elementary_keqs_iqr
from src.kinetics_integration.analyze_entropy import calculate_enzyme_entropy
from src.kinetics_integration.analyze_time_courses import time_course_analysis
from src.kinetics_integration.plot_clustermaps import cluster_map_param_inf
import os


def run_fitness_analysis(base_folder, file_in_base, model_type_list, ssd_threshold, Keq, column_order,
                         column_labels, x_label, scale_data, convert_to_ratios, n_model_ensembles):


    color_list = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']

    if not os.path.exists(''.join([base_folder, 'plots'])):
        os.makedirs(''.join([base_folder, 'plots']))

    Keq_local = None
    file_out_base = ''.join([base_folder, 'plots/TALA2_fitness_dist_all_', param_inf_type])
    plot_fitness_dist(file_in_base, file_out_base, model_type_list, n_model_ensembles, ssd_threshold, Keq_local,
                      column_order, column_labels, x_label, convert_to_ratios, color_list, filter=False)

    Keq_local = Keq
    file_out_base = ''.join([base_folder, 'plots/TALA2_n_valid_models_Keq_check_', param_inf_type])
    plot_n_good_models_dist(file_in_base, file_out_base, model_type_list, n_model_ensembles, ssd_threshold, Keq_local, column_order,
                            column_labels, x_label, scale_data, convert_to_ratios, color_list, filter=True)

    Keq_local = None
    file_out_base = ''.join([base_folder, 'plots/TALA2_n_valid_models_NO_Keq_check_', param_inf_type])
    plot_n_good_models_dist(file_in_base, file_out_base, model_type_list, n_model_ensembles, ssd_threshold, Keq_local, column_order,
                            column_labels, x_label, scale_data, convert_to_ratios, color_list, filter=True)


def analyze_talb(base_folder, enzyme, ssd_threshold, Keq, e_total, clustermaps, fitness_analysis, keq_range_analysis, entropy_analysis,
                 time_courses, param_inf_type):


    scale_data = False
    column_order = [0, 3, 5, 2, 4, 1]
    n_variables = len(column_order)

    convert_to_ratios = True

    model_type_list = ['all', 'dKd',  'Keq', 'Km1', 'Km2', 'Km3', 'Km4', 'kcat', 'Ki']
    column_labels = ['None', '$\Delta K_b$', '$K_{eq}$', '$K_m^{G3P}$', '$K_m^{E4P}$', '$K_m^{F6P}$', '$K_m^{S7P}$', '$k_{cat}$', '$K_i$']
    x_label = 'Data point removed'

    file_in_base = ''.join([base_folder, 'treated_data/rateconst_TALA2_'])

    if clustermaps:
        print 'plotting clustermaps'
        file_out_base = ''.join([base_folder, 'clustermaps/TALA2_ratios_', param_inf_type])

        if not os.path.exists(''.join([base_folder, 'clustermaps'])):
            os.makedirs(''.join([base_folder, 'clustermaps']))

        vmin = -6
        vmax = 9

        model_types = ['all', 'dKd']
        n_model_ensembles = 1

        Keq_local = None
        #cluster_map_param_inf(file_in_base, file_out_base, ssd_threshold, convert_to_ratios, column_order, vmin, vmax,
        #                      model_types, n_model_ensembles, Keq=Keq_local)

        fig_size = (8, 6)
        file_out_base = ''.join([base_folder, 'clustermaps/TALA2_rateconst_', param_inf_type])
        column_order_rateconst = [0, 1, 6, 7, 10, 11, 4, 5, 8, 9, 2, 3]
        cluster_map_param_inf(file_in_base, file_out_base, ssd_threshold, False, column_order_rateconst, vmin, vmax,
                              model_types, n_model_ensembles, fig_size, Keq=Keq_local)

    if fitness_analysis:
        print 'analyzing model fitness'
        n_model_ensembles = 100
        run_fitness_analysis(base_folder, file_in_base, model_type_list, ssd_threshold, Keq, column_order,
                             column_labels, x_label, scale_data, convert_to_ratios, n_model_ensembles)

    if keq_range_analysis:
        print 'analyzing var IQR'

        if param_inf_type == 'Keq':
            zorder_list = [7, 6, 0, 1, 2, 3, 4, 5]
            legend_n_cols = 4

        else:
            zorder_list = [8, 7, 0, 1, 2, 3, 4, 5, 6]
            legend_n_cols = 5

        file_out_base = ''.join([base_folder, 'plots/TALA2_elementary_rateconst_range_', param_inf_type])

        color_list = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']

        legend_details = [(0.5, 1.25), 0.82]


        x_labels = ['$\overleftarrow{k_{1}}$', '$\overrightarrow{k_{1}}$', '$\overleftarrow{k_{2}}$',
                    '$\overrightarrow{k_{2}}$', '$\overleftarrow{k_{3}}$', '$\overrightarrow{k_{3}}$',
                    '$\overleftarrow{k_{4}}$', '$\overrightarrow{k_{4}}$', '$\overleftarrow{k_{5}}$',
                    '$\overrightarrow{k_{5}}$', '$\overleftarrow{k_{6}}$', '$\overrightarrow{k_{6}}$']
        print len(x_labels)

        column_order_rateconst = [0, 1, 6, 7, 10, 11, 4, 5, 8, 9, 2, 3]
        n_model_ensembles = 100
        Keq_local = None
        calculate_elementary_keqs_iqr(file_in_base, file_out_base, model_type_list, n_model_ensembles, ssd_threshold,
                                      Keq_local, column_order_rateconst, False, column_labels, x_labels, color_list,
                                      zorder_list)
    if entropy_analysis:

        y_lims = [0.13, 0.19]

        if not os.path.exists(''.join([base_folder, 'entropy'])):
            os.makedirs(''.join([base_folder, 'entropy']))

        color_list = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']

        file_in_base = ''.join([base_folder, 'treated_data/rateconst_TALA2_'])

        plot_entropy_only = False
        n_model_ensembles = 100
        bin_width_list = [3]
        fixed_size = False
        Keq_local = None

        """
        column_order = [0, 3, 5, 2, 4, 1]

        file_out_base = ''.join([base_folder, 'entropy/TALA2_ratios_'])
        limit = 100
        calculate_enzyme_entropy(file_in_base, file_out_base, n_model_ensembles, model_type_list, column_order,
                                 n_variables, ssd_threshold, Keq_local, convert_to_ratios, limit, bin_width_list,
                                 column_labels, y_lims, x_label, color_list, plot_entropy_only=plot_entropy_only,
                                 fixed_size=fixed_size, n_samples_per_bin=None)
        limit = 60
        calculate_enzyme_entropy(file_in_base, file_out_base, n_model_ensembles, model_type_list, column_order,
                                 n_variables, ssd_threshold, Keq_local, convert_to_ratios, limit, bin_width_list,
                                 column_labels, y_lims, x_label, color_list, plot_entropy_only=plot_entropy_only,
                                 fixed_size=fixed_size, n_samples_per_bin=None)
        """

        file_out_base = ''.join([base_folder, 'entropy/TALA2_rateconst_'])
        column_order_rateconst = [0, 1, 6, 7, 10, 11, 4, 5, 8, 9, 2, 3]
        n_variables_rateconst = len(column_order_rateconst)
        limit = 100
        calculate_enzyme_entropy(file_in_base, file_out_base, n_model_ensembles, model_type_list,
                                 column_order_rateconst, n_variables_rateconst, ssd_threshold, Keq_local, False, limit,
                                 bin_width_list, column_labels, y_lims, x_label, color_list,
                                 plot_entropy_only=plot_entropy_only, fixed_size=fixed_size, n_samples_per_bin=None)
        limit = 30
        calculate_enzyme_entropy(file_in_base, file_out_base, n_model_ensembles, model_type_list,
                                 column_order_rateconst, n_variables_rateconst, ssd_threshold, Keq_local, False, limit,
                                 bin_width_list, column_labels, y_lims, x_label, color_list,
                                 plot_entropy_only=plot_entropy_only, fixed_size=fixed_size, n_samples_per_bin=None)

    if time_courses:

        n_model_ensembles = 100
        plot_only = False

        subs_list = ['s7p' 'g3p']
        prod_list = ['f6p', 'e4p']
        perturbation_list = ['noPerturb']

        filter_models = False
        x_lims = [0, 100]

        time_point_list = [0, 100]
        plot_type_list = ['mets'] #'enz', 'mets', 'flux']

        for perturb in perturbation_list:

            for plot_type in plot_type_list:

                mets_order = None
                legend_label = None

                if plot_type == 'mets':

                    file_type = 'conc_mets'
                    #y_lims = [0, 0.0025]  # 0.00022
                    species_list = ['e4p', 'f6p', 'g3p', 's7p']
                    mets_order = ['f6p', 'e4p', 'g3p', 's7p']
                    y_label = 'Concentration (mol/L)'
                    y_lims = None #[[0.0000030, 0.0000045], [0.000220, 0.000235], [0.0014, 0.0016], [0.000085, 0.000092], [0.00232, 0.00235]]
                    #y_lims_timecourses = [[5.95*10**-4, 6.15*10**-4], [9.5*10**-5, 1.15*10**-4], [2.05*10**-4, 2.2*10**-4],
                    #                      [2.6*10**-4, 2.8*10**-4], [2.32*10**-3, 2.34*10**-3]]
                    y_lims_timecourses = [[2.4*10**-3, 2.6*10**-3], [5*10**-5, 10**-4], [1.9*10**-4, 3*10**-4],
                                           [8*10**-4, 9*10**-4]]
                    #y_lims_timecourses = None
                    plot_spacings = [0.89, 0.15, 0.16, 0.98]
                    fig_size = (6,4)

                elif plot_type == 'enz':
                    file_type = 'conc_enz'
                    #y_lims = [0, 7 * 10 ** -5]
                    species_list = ['TALA2', 'TALA2&f6p', 'TALA2&mod', 'TALA2&pi', 'TALA2&s7p', 'TALA2&mod&e4p', 'TALA2&mod&g3p', 'TALA2&mod&pi']
                    #species_list = ['GAPD', 'GAPD&nad', 'GAPD&nadh', 'GAPD&nad&g3p', 'GAPD&nadh&13dpg', 'GAPD&nad&g3p&pi']
                    legend_label = ['E', 'E-f6p', 'E\'',  'E\'-e4p', 'E\'-g3p', 'E\'-pi', 'E-pi', 'E-s7p']
                    y_label = 'Concentration (mol/L)'
                    y_lims = None # [[0, 0.00004], [0, 0.000008],[0, 0.00008], [0., 0.00008], [0., 0.00003], [0., 0.00003]]
                    y_lims_timecourses = None
                    plot_spacings = [0.89, 0.2, 0.15, 0.98]
                    fig_size = (6,4)

                elif plot_type == 'flux':
                    file_type = 'flux'
                    #y_lims = [-0.001, 0.001]
                    species_list = ['v_tala21', 'v_tala22', 'v_tala23', 'v_tala24', 'v_tala25', 'v_tala26']
                    legend_label = ['tala21', 'tala22', 'tala23', 'tala24', 'tala25', 'tala26']
                    y_label = 'Flux (mol/L/s)'
                    y_lims = None
                    y_lims_timecourses = None
                    plot_spacings = [0.89, 0.15, 0.15, 0.98]
                    fig_size = (6,4)

                else:
                    print 'plot_type not specified correctly, ', plot_type
                    return

                file_in_base = ''.join([base_folder, 'model_simulations/data/', perturb, '/sim_res_', file_type])
                file_out_base = ''.join([base_folder, 'model_simulations/plots/', perturb, '/', plot_type, '/sim_res_', file_type])

                if not os.path.exists(''.join([base_folder, 'model_simulations/plots/', perturb, '/', plot_type])):
                    os.makedirs(''.join([base_folder, 'model_simulations/plots/', perturb, '/', plot_type]))

                time_course_analysis(file_in_base, file_out_base, species_list, Keq, subs_list, prod_list, mets_order,
                                     legend_label, n_model_ensembles, plot_type, time_point_list, plot_spacings,
                                     fig_size, y_lims=y_lims, y_lims_timecourses=y_lims_timecourses,
                                     plot_only=plot_only)


if __name__ == '__main__' :
    enzyme = 'TALA2'
    ssd_threshold = 0.1
    Keq = 1.31
    #Keq = ''
    param_inf_type = ''
    e_total = 5.4*10**-6

    clustermaps = False
    fitness_analysis = True
    keq_range_analysis = False
    entropy_analysis = False
    time_courses = False


    base_folder = '/home/mrama/Desktop/MD/eMASS-MD_complete_data/enzyme_models/TALA2/TALA2_param_inf2/output/'

    analyze_talb(base_folder, enzyme, ssd_threshold, Keq, e_total, clustermaps, fitness_analysis, keq_range_analysis,
                 entropy_analysis, time_courses, param_inf_type)