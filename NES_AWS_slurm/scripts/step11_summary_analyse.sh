#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
    echo "\$PROJ_HOME/config.bash not found"
    exit
fi
source ${PROJ_HOME}config.bash

source /share/local/miniforge3/etc/profile.d/conda.sh
conda activate AZtutorial
source ${PROJ_HOME}config.bash

results_csv="${PROJ_HOME}workpath/summary/results.csv"

python ${SCRIPTS_PATH}results_analyse.py ${PROJ_HOME} ${results_csv}
