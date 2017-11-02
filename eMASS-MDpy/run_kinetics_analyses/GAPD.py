from src.kinetics_integration.analyze_model_fitness import plot_fitness_dist, plot_n_good_models_dist
from src.kinetics_integration.analyze_var_IQR import calculate_elementary_keqs_iqr
from src.kinetics_integration.analyze_entropy import calculate_enzyme_entropy
from src.kinetics_integration.analyze_time_courses import time_course_analysis
from src.kinetics_integration.plot_clustermaps import cluster_map_param_inf
import os


def run_fitness_analysis(base_folder, file_in_base, model_type_list, ssd_threshold, Keq, column_order,
                         column_labels, x_label, scale_data, convert_to_ratios, n_model_ensembles, param_inf_type):

    color_list = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']

    if not os.path.exists(''.join([base_folder, 'plots'])):
        os.makedirs(''.join([base_folder, 'plots']))

    Keq_local = None
    file_out_base = ''.join([base_folder, 'plots/GAPD_fitness_dist_all_', param_inf_type])
    plot_fitness_dist(file_in_base, file_out_base, model_type_list, n_model_ensembles, ssd_threshold, Keq_local,
                      column_order, column_labels, x_label, convert_to_ratios, color_list, filter=False)

    Keq_local = Keq
    file_out_base = ''.join([base_folder, 'plots/GAPD_valid_models_Keq_check_', param_inf_type])
    plot_n_good_models_dist(file_in_base, file_out_base, model_type_list, n_model_ensembles, ssd_threshold, Keq_local,
                            column_order, column_labels, x_label, scale_data, convert_to_ratios, color_list, filter=True)


    Keq_local = None
    file_out_base = ''.join([base_folder, 'plots/GAPD_n_valid_models_noKeq_NO_Keq_check_', param_inf_type])
    plot_n_good_models_dist(file_in_base, file_out_base, model_type_list, n_model_ensembles, ssd_threshold, Keq_local,
                            column_order, column_labels, x_label, scale_data, convert_to_ratios, color_list, filter=True)


def analyze_gapd(base_folder, enzyme, Keq, e_total, clustermaps, fitness_analysis, keq_range_analysis, entropy_analysis,
                 time_courses, param_inf_type):

    ssd_threshold = 1
    scale_data = False

    column_order = [0, 2, 4, 5, 3, 1]
    n_variables = len(column_order)
    convert_to_ratios = True

    model_type_list = ['all',  'dKd', 'Keq', 'Kd', 'Km1', 'Km2', 'Km3', 'kcat']
    column_labels = ['None', '$\Delta K_d$', '$K_{eq}$', '$K_d^{nad}$', '$K_m^{nad}$', '$K_m^{g3p}$', '$K_m^{pi}$', '$k_{cat}$']
    x_label = 'Data point removed'

    file_in_base = ''.join([base_folder, 'treated_data/rateconst_GAPD_'])

    if clustermaps:
        print 'plotting clustermaps'

        file_out_base = ''.join([base_folder, 'clustermaps/GAPD_ratio_', param_inf_type])

        if not os.path.exists(''.join([base_folder, 'clustermaps'])):
            os.makedirs(''.join([base_folder, 'clustermaps']))

        vmin = -6
        vmax = 9

        model_types = ['all', 'dKd']
        n_model_ensembles = 100

        Keq_local = None
        #cluster_map_param_inf(file_in_base, file_out_base, ssd_threshold, convert_to_ratios, column_order, vmin, vmax,
        #                      model_types, n_model_ensembles, Keq=Keq_local)


        fig_size = (8, 6)
        file_out_base = ''.join([base_folder, 'clustermaps/GAPD_rateconst_', param_inf_type])
        column_order_rateconst = [0, 1, 4, 5, 8, 9, 10, 11, 6, 7, 2, 3]
        cluster_map_param_inf(file_in_base, file_out_base, ssd_threshold, False, column_order_rateconst, vmin, vmax,
                              model_types, n_model_ensembles, fig_size, Keq=Keq_local)

    if fitness_analysis:
        print 'analyzing model fitness'

        n_model_ensembles = 100
        run_fitness_analysis(base_folder, file_in_base, model_type_list, ssd_threshold, Keq, column_order,
                             column_labels, x_label, scale_data, convert_to_ratios, n_model_ensembles, param_inf_type)

    if keq_range_analysis:
        print 'analyzing var IQR'

        if param_inf_type == 'Keq':
            zorder_list = [6, 0, 1, 2, 3, 4, 5]
        else:
           zorder_list = [7, 6, 0, 1, 2, 3, 4, 5]
        print zorder_list

        file_out_base = ''.join([base_folder, 'plots/GAPD_elementary_rateconst_range_', param_inf_type])

        x_labels = ['$\overleftarrow{k_{1}}$', '$\overrightarrow{k_{1}}$', '$\overleftarrow{k_{2}}$',
                    '$\overrightarrow{k_{2}}$', '$\overleftarrow{k_{3}}$', '$\overrightarrow{k_{3}}$',
                    '$\overleftarrow{k_{4}}$', '$\overrightarrow{k_{4}}$', '$\overleftarrow{k_{5}}$',
                    '$\overrightarrow{k_{5}}$', '$\overleftarrow{k_{6}}$', '$\overrightarrow{k_{6}}$']

        color_list = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']

        legend_details = [(0.5, 1.25), 0.82]
        legend_n_cols = 4
        column_order_rateconst = [0, 1, 4, 5, 8, 9, 10, 11, 6, 7, 2, 3]

        Keq_local = None
        n_model_ensembles = 100

        calculate_elementary_keqs_iqr(file_in_base, file_out_base, model_type_list, n_model_ensembles, ssd_threshold,
                                      Keq_local, column_order_rateconst, False, column_labels, x_labels, color_list,
                                      zorder_list)

    if entropy_analysis:
        print 'analyzing entropy'

        y_lims = [0.15, 0.25]

        if not os.path.exists(''.join([base_folder, 'entropy'])):
            os.makedirs(''.join([base_folder, 'entropy']))

        # x_labels = ['', '$K_1$', '$K_2$', '$K_3$',  '$K_4$', '$K_5$', '$K_6$', '']
        color_list = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']

        plot_entropy_only = False
        n_model_ensembles = 100
        bin_width_list = [3]
        fixed_size = False
        Keq_local = None
        limit = 100
        """
        file_out_base = ''.join([base_folder, 'entropy/GAPD_ratios_'])
        calculate_enzyme_entropy(file_in_base, file_out_base, n_model_ensembles, model_type_list, column_order,
                                 n_variables, ssd_threshold, Keq_local, convert_to_ratios, limit, bin_width_list,
                                 column_labels, y_lims, x_label, color_list, plot_entropy_only=plot_entropy_only,
                                 fixed_size=fixed_size, n_samples_per_bin=None)
        """

        file_out_base = ''.join([base_folder, 'entropy/GAPD_rateconst_'])
        column_order_rateconst = [0, 1, 4, 5, 8, 9, 10, 11, 6, 7, 2, 3]
        n_variables_rateconst = len(column_order_rateconst)
        calculate_enzyme_entropy(file_in_base, file_out_base, n_model_ensembles, model_type_list,
                                 column_order_rateconst, n_variables_rateconst, ssd_threshold, Keq_local, False, limit,
                                 bin_width_list, column_labels, y_lims, x_label, color_list,
                                 plot_entropy_only=plot_entropy_only, fixed_size=fixed_size, n_samples_per_bin=None)

    if time_courses:
        print 'analyzing time-courses'

        n_model_ensembles = 100
        plot_only = False

        subs_list = ['g3p', 'nad', 'pi']
        prod_list = ['13dpg', 'nadh']
        perturbation_list = ['noPerturb']

        filter_models = False
        x_lims = [0, 10]

        time_point_list = [0, 10]
        plot_type_list = ['enz', 'mets', 'flux']

        for perturb in perturbation_list:
            for plot_type in plot_type_list:
                mets_order = None
                legend_label = None

                if plot_type == 'mets':

                    file_type = 'conc_mets'
                    #y_lims = [0, 0.0025]  # 0.00022
                    species_list = ['13dpg', 'g3p', 'nad', 'nadh', 'pi']
                    mets_order = ['g3p', 'nad', 'pi', '13dpg', 'nadh']
                    y_label = 'Concentration (mol/L)'
                    y_lims = [[0.0000030, 0.0000045], [0.000220, 0.000235], [0.0014, 0.0016], [0.000085, 0.000092], [0.00232, 0.00235]]
                    y_lims_timecourses = [[2.42*10**-3, 2.5*10**-3], [2.38*10**-2, 2.39*10**-2],
                                          [7*10**-5, 1.22*10**-4], [1.8*10**-4, 2.5*10**-4], [4*10**-5, 8*10**-5]]
                    # y_lims_timecourses = None
                    y_lims=None
                    plot_spacings = [0.89, 0.15, 0.15, 0.98]
                    fig_size = (6,4)

                elif plot_type == 'enz':
                    file_type = 'conc_enz'
                    #y_lims = [0, 7 * 10 ** -5]
                    #species_list = ['GAPD', 'GAPD&13dpg', 'GAPD&nad', 'GAPD&13dpg&nadh', 'GAPD&nad&g3p', 'GAPD&nad&g3p&pi']
                    species_list = ['GAPD', 'GAPD&nad', 'GAPD&nadh', 'GAPD&nad&g3p', 'GAPD&nadh&13dpg', 'GAPD&nad&g3p&pi']
                    legend_label = ['E', 'E-nad',  'E-nad-g3p', 'E-nad-g3p-pi', 'E-nadh', 'E-nadh-13dpg']
                    y_label = 'Concentration (mol/L)'
                    y_lims = [[0, 0.00004], [0, 0.000008],[0, 0.00008], [0., 0.00008], [0., 0.00003], [0., 0.00003]]
                    y_lims_timecourses = None
                    plot_spacings = [0.89, 0.3, 0.15, 0.98]
                    fig_size=(6,5)

                elif plot_type == 'flux':
                    file_type = 'flux'
                    #y_lims = [-0.001, 0.001]
                    species_list = ['v_gapd1', 'v_gapd2', 'v_gapd3', 'v_gapd4', 'v_gapd5', 'v_gapd6']
                    legend_label = ['gapd1', 'gapd2', 'gapd3', 'gapd4', 'gapd5', 'gapd6']
                    y_label = 'Flux (mol/L/s)'
                    y_lims = None
                    y_lims_timecourses = None
                    plot_spacings = [0.89, 0.15, 0.2, 1]
                    fig_size= (6,4)

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
    enzyme = 'GAPD'
    Keq = 0.452  # to filter by keq
    #Keq = ''
    param_inf_type = ''
    e_total = 6.875*10**-5

    clustermaps = True
    fitness_analysis = True
    keq_range_analysis = True
    entropy_analysis = True
    time_courses = False

    base_folder = '/home/mrama/Dropbox/PhD_stuff/Projects/MD/eMASS-MD/enzyme_models/GAPD/GAPD_param_inf/output/'

    analyze_gapd(base_folder, enzyme, Keq, e_total, clustermaps, fitness_analysis, keq_range_analysis, entropy_analysis,
                 time_courses, param_inf_type)