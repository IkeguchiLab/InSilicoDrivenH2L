#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
    echo "\$PROJ_HOME/config.bash not found"
    exit
fi
source ${PROJ_HOME}config.bash

source /share/local/miniforge3/etc/profile.d/conda.sh

conda activate AZtutorial
source ${PROJ_HOME}config.bash

transition_dir=`find ${TRANSITION_PATH}edge*/stateA/run*_* -type d -name "run*_*"`
transition_dir_array=(${transition_dir})
len_transition_dir_array=${#transition_dir_array[@]}

for ((i=0; i<len_transition_dir_array; i++));do
dir_temp=${transition_dir_array[i]}
run_name=$(basename ${dir_temp})
echo ${run_name}
IFS='_' read -r -a parts <<< "${run_name}"
run_number="${parts[0]}"
protein_water="${parts[1]}"
echo ${run_number}
echo ${protein_water}

cd ${transition_dir_array[i]}
mkdir -p ../../analyse_${run_name}

analysis_dir="../../analyse_${run_name}"
results_file="${analysis_dir}/results.txt"

# Check if results.txt already exists
if [ -f "${results_file}" ]; then
    echo "Results file already exists for ${run_name}, skipping..."
    continue
fi

pmx analyse -fA *xvg -fB ../../stateB/${run_name}/*xvg -o ../../analyse_${run_name}/result.txt -oA ../../analyse_${run_name}/integ0.dat -oB ../../analyse_${run_name}/integ1.dat -w ../../analyse_${run_name}/wplot.png -o ../../analyse_${run_name}/results.txt -t 298 -b 100 --units kcal --reverseB

done
