#!/usr/bin/env bash
# this script should be run from the R subdirectory

config_file="../../experiments/all_configurations.csv"
num_configs=`wc -l < $config_file`
echo $num_configs
for i in `seq 1 $num_configs`
do
    cmd="qsub -cwd -v PATH,R_HOME run_instance.R $i $config_file ../../experiments/results.csv"
    echo $cmd
    bash -c "$cmd"
done

