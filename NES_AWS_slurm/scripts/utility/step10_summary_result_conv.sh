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

python ${SCRIPTS_PATH}results_txt_conv.py ${PROJ_HOME}
