#!/usr/bin/env bash
# this script should be run from the R subdirectory

config_file="../../experiments/all_configurations_rw.csv"
num_configs=`wc -l < $config_file`
echo $num_configs
for i in `seq 1 $num_configs`
do
    cmd="qsub -cwd -o logs/out${i}.txt -e logs/error${i}.txt -l mem_free=8G -l mem_token=8G -v PATH,R_HOME run_instance.R $i $config_file ../../experiments/results_rw.csv"
    echo $cmd
    bash -c "$cmd"
done


