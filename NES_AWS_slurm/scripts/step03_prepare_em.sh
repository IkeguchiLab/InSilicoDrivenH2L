#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
    echo "\$PROJ_HOME/config.bash not found"
    exit
fi
source ${PROJ_HOME}config.bash

while read lig
do
run_dir=`find ${EQ_PATH}${lig}/ -type d -name "run*" | sort`
run_dir_array=(${run_dir})
for run_directory in "${run_dir_array[@]}";do
echo ${run_directory}
mkdir -p ${run_directory}/cal01_em_protein
${GMX} grompp -f ${MDP_PATH}em_l0.mdp -c ${EQ_PATH}${lig}/protein/ions.pdb -r ${EQ_PATH}${lig}/protein/ions.pdb -p ${EQ_PATH}${lig}/protein/topol.top -o ${run_directory}/cal01_em_protein/tpr.tpr -maxwarn 4 -po ${run_directory}/cal01_em_protein/mdout.mdp >& ${run_directory}/cal01_em_protein/log_prepare_em
sleep 0.2
mkdir -p ${run_directory}/cal01_em_water
${GMX} grompp -f ${MDP_PATH}em_l0.mdp -c ${EQ_PATH}${lig}/water/ions.pdb -r ${EQ_PATH}${lig}/water/ions.pdb -p ${EQ_PATH}${lig}/water/topol.top -o ${run_directory}/cal01_em_water/tpr.tpr -maxwarn 4 -po ${run_directory}/cal01_em_water/mdout.mdp >& ${run_directory}/cal01_em_water/log_prepare_em
done
done < "${LIGAND_PATH}liglist.txt"
