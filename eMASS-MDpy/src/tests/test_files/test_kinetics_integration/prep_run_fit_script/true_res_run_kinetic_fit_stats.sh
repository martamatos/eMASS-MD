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

start_fit=1
end_fit=$((100-1))

for i in `seq $start_fit $end_fit`;
do
	echo "Fit Number: $i"
	python $scripts_dir/run_fit.py $data_dir/psoParameters.txt $data_dir/lmaParameters.txt $results_dir/summary.txt $results_dir/psoResults_${i}.txt $results_dir/lmaResults_${i}.txt $num_trials $data_dir/G6PDH2r.dat
done
