#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
    echo "\$PROJ_HOME/config.bash not found"
    exit
fi
source ${PROJ_HOME}config.bash

results_file=`find ${PROJ_HOME}workpath/transition/edge_*/analyse_*/results.txt -type f`
results_file_array=(${results_file})
len_results_file_array=${#results_file_array[@]}

mkdir -p ${PROJ_HOME}workpath/summary

#for ((i=0; i<len_results_file_array; i++));do
#echo ${results_file_array[i]}
#python ${SCRIPTS_PATH}results_txt.py ${results_file_array[i]}
#done

python ${SCRIPTS_PATH}results_txt.py ${PROJ_HOME}
