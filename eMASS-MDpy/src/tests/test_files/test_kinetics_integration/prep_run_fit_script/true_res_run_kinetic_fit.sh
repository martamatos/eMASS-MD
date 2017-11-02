#!/bin/sh
#PBS -N G6PDH2r2
#PBS -q hpc
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=1,mem=2gb

export PYTHONPATH=/zhome/89/0/90554/.local/lib/python2.6/site-packages

main_dir=/zhome/89/0/90554/kinetics_fit/G6PD
data_dir=/zhome/89/0/90554/kinetics_fit/G6PD/data
results_dir=/zhome/89/0/90554/kinetics_fit/G6PD/results
scripts_dir=/zhome/89/0/90554/kinetics_fit/scripts

header="best.fitness, num_generations, pop_size, neighborhood_size, inertia, cognitive_rate, social_rate"
echo $header > $results_dir/summary_$PBS_ARRAYID.txt

num_trials=$((20-1))

datasetList=(Km Keq KmKeqdKd1 KmKeqdKd2)
n_datasets=$((${#datasetList[*]}-1))
echo $datasetList

dKd_list=(1.e-9 1.e-8 1.e-7 0.000001 0.00001 0.0001 0.001 0.01 0.1 1. 10 100 1000 10000 100000 1000000 10000000 100000000 1000000000)
n_dKd=$((${#dKd_list[*]}-1))
echo $dKd_list

for i in `seq 0 $n_datasets`;
do
	echo $i
	echo ${datasetList[$i]}

	for j in `seq 0 $n_dKd`;
	do
		echo $j
		echo ${dKd_list[$j]}
		echo $main_dir/G6PDH2r_${datasetList[$i]}_${dKd_list[$j]}.dat

		python $scripts_dir/run_fit.py $data_dir/psoParameters.txt $data_dir/lmaParameters.txt $results_dir/summary_${datasetList[$i]}_${dKd_list[$j]}.txt $results_dir/psoResults_${datasetList[$i]}_${dKd_list[$j]}.txt $results_dir/lmaResults_${datasetList[$i]}_${dKd_list[$j]}.txt $num_trials $data_dir/G6PDH2r_${datasetList[$i]}_${dKd_list[$j]}.dat
	done
done
