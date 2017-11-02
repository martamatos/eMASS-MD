import os
import subprocess
import time


def _write_header(f_out, main_dir, data_file, enzyme_folder, python_packages_path, job_number, n_trials, walltime):

    f_out.write('#!/bin/sh\n')
    f_out.write(''.join(['#PBS -N ', data_file, str(job_number), '\n']))
    f_out.write('#PBS -q hpc\n') # user will probably need to change the queue here
    f_out.write(''.join(['#PBS -l walltime=', walltime, '\n']))
    f_out.write('#PBS -l nodes=1:ppn=1,mem=2gb\n\n')

    f_out.write(''.join(['export PYTHONPATH=', python_packages_path, '\n\n']))

    f_out.write(''.join(['main_dir=', main_dir, enzyme_folder, '\n']))
    f_out.write(''.join(['data_dir=', main_dir, enzyme_folder, '/data\n']))
    f_out.write(''.join(['results_dir=', main_dir, enzyme_folder, '/results\n']))
    f_out.write(''.join(['scripts_dir=', main_dir, 'scripts\n\n']))

    f_out.write(
        'header="best.fitness, num_generations, pop_size, neighborhood_size, inertia, cognitive_rate, social_rate"\n')
    f_out.write('echo $header > $results_dir/summary_$PBS_ARRAYID.txt\n\n')

    f_out.write(''.join(['num_trials=', str(n_trials), '\n\n']))


def prep_run_fit_script(file_out, main_dir, data_file, enzyme_folder, python_packages_path, job_number, n_trials,
                        start_fit, end_fit, walltime='5:00:00', qsub=False, fit_label=''):
    """
        Generates a bash script to run the same enzyme fit (end_fit - start_fit) times.
        This bash script is basically a cluster job to be executed on a cluster.

    :param file_out: name of the bash script to be generated to run the fits.
    :param main_dir: main folder for the fit results.
    :param data_file: file name with simulated data.
    :param enzyme_folder: name of the folder where the enzyme data is and the fit results will be stored.
    :param python_packages_path: path to where you have
    :param job_number: number of the job to be submitted.
    :param n_trials: how many times the fit is repeated.
    :param start_fit: the number of the first fit
    :param end_fit: the number of the last fit
    :param walltime: walltime for the job.
    :param qsub: a boolean specifying whether or not the job should be submitted automatically.
    :param fit_label: label for the fit being done, to be attached to every result file name.
    :return: None
    """

    with open(file_out, 'w') as f_out:
        _write_header(f_out, main_dir, data_file, enzyme_folder, python_packages_path, job_number, n_trials, walltime)

        f_out.write(''.join(['fit_label=', str(fit_label), '\n']))

        f_out.write(''.join(['start_fit=', str(start_fit), '\n']))
        f_out.write(''.join(['end_fit=', str(end_fit), '\n\n']))
        f_out.write('start_t=$SECONDS\n\n')

        f_out.write('for i in `seq $start_fit $end_fit`;\n')
        f_out.write('do\n')
        f_out.write('\techo "Fit Number: $i"\n')

        f_out.write(''.join(['\tpython $scripts_dir/run_fit_rel.py $data_dir/psoParameters.txt ',
                             '$data_dir/lmaParameters.txt $results_dir/summary.txt '
                             '$results_dir/psoResults_${fit_label}_${i}.txt $results_dir/lmaResults_${fit_label}_${i}.txt ',
                             '$num_trials $data_dir/', data_file, '.dat\n']))

        f_out.write('done\n\n')
        f_out.write('end_t=$SECONDS\n')
        f_out.write('echo "time elapsed $(($end_t - $start_t))"\n')

    os.chmod(file_out, 0755)
    if qsub:
        subprocess.call(''.join(['qsub ', file_out]), shell=True)


def run_ENO_param_scan_ensemble():
    # path to folder where all data and scripts are stored
    main_dir = '/zhome/89/0/90554/kinetics_fit/'
    # path to folder where python packages are installed, important if they are only installed in the user's folder
    python_packages_path = '/zhome/89/0/90554/.local/lib/python2.6/site-packages'
    enzyme = 'ENO'
    enzyme_folder = 'ENO_param_scan'
    dataset = 'customRatio_1'

    dKd_list = ['0.000000000001', '0.00000000001',
                '0.0000000001', '0.000000001',
                '0.00000001', '0.0000001',
                '0.000001', '0.00001', '0.0001',
                '0.001', '0.01', '0.1',
                '0.00000175', '1.', '10.',
                '100.', '1000.', '10000.',
                '100000.', '1000000.', '10000000.',
                '100000000.', '1000000000.', '10000000000.',
                '100000000000.', '1000000000000.']

    n_trials = 100
    walltime = '24:00:00'
    qsub = True
    n_ensembles = 1

    for i, dKd in enumerate(dKd_list):
        fit_label = '_'.join([dataset, dKd])
        data_file = '_'.join([enzyme, fit_label])
        job_number = i
        start_fit = 1
        end_fit = n_ensembles
        file_out = ''.join(['/zhome/89/0/90554/run_kinetic_fit_', data_file, '_', str(job_number), '.sh'])
        prep_run_fit_script(file_out, main_dir, data_file, enzyme_folder, python_packages_path, job_number,
                            n_trials, start_fit, end_fit, walltime=walltime, qsub=qsub, fit_label=fit_label)
        time.sleep(2)


def run_ENO_param_influence():
    # path to folder where all data and scripts are stored
    main_dir = '/zhome/89/0/90554/kinetics_fit/'
    # path to folder where python packages are installed, important if they are only installed in the user's folder
    python_packages_path = '/zhome/89/0/90554/.local/lib/python2.6/site-packages'
    enzyme = 'ENO'
    enzyme_folder = 'ENO_param_inf'
    dataset_list = ['all', 'Keq', 'dKd', 'Km', 'kcat']

    n_trials = 100
    walltime = '24:00:00'
    qsub = True

    n_ensembles = 100
    n_models_block = 10
    n_jobs = n_ensembles / n_models_block

    for fit_label in dataset_list:
        data_file = '_'.join([enzyme, fit_label])
        for i in range(n_jobs):
            job_number = i
            start_fit = i * n_models_block
            end_fit = (i + 1) * n_models_block - 1
            file_out = ''.join(['/zhome/89/0/90554/run_kinetic_fit_', data_file, '_', str(job_number), '.sh'])
            prep_run_fit_script(file_out, main_dir, data_file, enzyme_folder, python_packages_path, job_number,
                                n_trials, start_fit, end_fit, walltime=walltime, qsub=qsub, fit_label=fit_label)
            time.sleep(2)


def run_GAPD_param_scan_ensemble():
    # path to folder where all data and scripts are stored
    main_dir = '/zhome/89/0/90554/kinetics_fit/'
    # path to folder where python packages are installed, important if they are only installed in the user's folder
    python_packages_path = '/zhome/89/0/90554/.local/lib/python2.6/site-packages'
    enzyme = 'GAPD'
    enzyme_folder = 'GAPD_param_scan'
    dataset = 'customRatio_1'

    dKd_list = ['0.000000000001', '0.00000000001',
                '0.0000000001', '0.000000001',
                '0.00000001', '0.0000001',
                '0.000001', '0.00001', '0.0001',
                '0.001', '0.01', '0.1',
                '0.00000109', '1.', '10.',
                '100.', '1000.', '10000.',
                '100000.', '1000000.', '10000000.',
                '100000000.', '1000000000.', '10000000000.',
                '100000000000.', '1000000000000.']


    n_trials = 100
    walltime = '24:00:00'
    qsub = True
    n_ensembles = 10

    for i, dKd in enumerate(dKd_list):
        fit_label = '_'.join([dataset, dKd])
        data_file = '_'.join([enzyme, fit_label])
        job_number = i
        start_fit = 1
        end_fit = n_ensembles
        file_out = ''.join(['/zhome/89/0/90554/run_kinetic_fit_', data_file, '_', str(job_number), '.sh'])
        prep_run_fit_script(file_out, main_dir, data_file, enzyme_folder, python_packages_path, job_number,
                            n_trials, start_fit, end_fit, walltime=walltime, qsub=qsub, fit_label=fit_label)
        time.sleep(2)


def run_GAPD_param_influence():
    # path to folder where all data and scripts are stored
    main_dir = '/zhome/89/0/90554/kinetics_fit/'
    # path to folder where python packages are installed, important if they are only installed in the user's folder
    python_packages_path = '/zhome/89/0/90554/.local/lib/python2.6/site-packages'
    enzyme = 'GAPD'
    enzyme_folder = 'GAPD_param_inf'
    dataset_list = ['all', 'Keq', 'Km1', 'Km2', 'Km3', 'kcat', 'Kd', 'dKd']

    n_trials = 100
    walltime = '48:00:00'
    qsub = True
    n_ensembles = 100
    n_models_block = 10
    n_jobs = n_ensembles / n_models_block

    for fit_label in dataset_list:
        data_file = '_'.join([enzyme, fit_label])
        for i in range(n_jobs):
            job_number = i
            start_fit = i * n_models_block
            end_fit = (i + 1) * n_models_block - 1
            file_out = ''.join(['/zhome/89/0/90554/run_kinetic_fit_', data_file, '_', str(job_number), '.sh'])
            prep_run_fit_script(file_out, main_dir, data_file, enzyme_folder, python_packages_path, job_number,
                                n_trials, start_fit, end_fit, walltime=walltime, qsub=qsub, fit_label=fit_label)
            time.sleep(2)


def run_TALA2_param_scan_ensemble():
    # path to folder where all data and scripts are stored
    main_dir = '/zhome/89/0/90554/kinetics_fit/'
    # path to folder where python packages are installed, important if they are only installed in the user's folder
    python_packages_path = '/zhome/89/0/90554/.local/lib/python2.6/site-packages'
    enzyme = 'TALA2'
    enzyme_folder = 'TALA2_param_scan'

    dataset = 'customRatio_1'

    dKd_list = ['0.000000000001', '0.00000000001',
                '0.0000000001', '0.000000001',
                '0.00000001', '0.0000001',
                '0.000001', '0.00001', '0.0001',
                '0.001', '0.01', '0.1',
                '1.', '10.', '7.33',
                '100.', '1000.', '10000.',
                '100000.', '1000000.', '10000000.',
                '100000000.', '1000000000.', '10000000000.',
                '100000000000.', '1000000000000.']

    n_trials = 100
    walltime = '24:00:00'
    qsub = True
    n_ensembles = 10

    for i, dKd in enumerate(dKd_list):
        fit_label = '_'.join([dataset, dKd])
        data_file = '_'.join([enzyme, fit_label])
        job_number = i
        start_fit = 1
        end_fit = n_ensembles
        file_out = ''.join(['/zhome/89/0/90554/run_kinetic_fit_', data_file, '_', str(job_number), '.sh'])
        prep_run_fit_script(file_out, main_dir, data_file, enzyme_folder, python_packages_path, job_number,
                            n_trials, start_fit, end_fit, walltime=walltime, qsub=qsub, fit_label=fit_label)
        time.sleep(2)


def run_TALA2_param_influence():
    # path to folder where all data and scripts are stored
    main_dir = '/zhome/89/0/90554/kinetics_fit/'
    # path to folder where python packages are installed, important if they are only installed in the user's folder
    python_packages_path = '/zhome/89/0/90554/.local/lib/python2.6/site-packages'
    enzyme = 'TALA2'
    enzyme_folder = 'TALA2_param_inf'
    dataset_list = ['all', 'dKd', 'Km1', 'Km2', 'Km3', 'Km4', 'kcat', 'Ki', 'Keq']


    n_trials = 100
    walltime = '12:00:00'
    qsub = True
    n_ensembles = 100
    n_models_block = 2
    n_jobs = n_ensembles / n_models_block

    for fit_label in dataset_list:
        data_file = '_'.join([enzyme, fit_label])
        for i in range(n_jobs):
            job_number = i
            start_fit = i * n_models_block
            end_fit = (i + 1) * n_models_block - 1
            file_out = ''.join(['/zhome/89/0/90554/run_kinetic_fit_', data_file, '_', str(job_number), '.sh'])
            prep_run_fit_script(file_out, main_dir, data_file, enzyme_folder, python_packages_path, job_number,
                                n_trials, start_fit, end_fit, walltime=walltime, qsub=qsub, fit_label=fit_label)
            time.sleep(2)


if __name__ == '__main__':
    run_ENO_param_scan_ensemble()
    #run_ENO_param_influence()
    #run_GAPD_param_scan_ensemble()
    #run_GAPD_param_influence()
    #run_TALA2_param_scan_ensemble()
    #run_TALA2_param_influence()
