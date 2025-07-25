#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
    echo "\$PROJ_HOME/config.bash not found"
    exit
fi
source ${PROJ_HOME}config.bash

edge_file=${TRANSITION_PATH}edgelist.txt

if [ -e "${edge_file}" ]; then
while IFS=',' read -r ligA ligB; do
#echo ${ligA}
#echo ${ligB}
edge_dir=edge_${ligA}_${ligB}
echo ${edge_dir}
mkdir -p ${TRANSITION_PATH}${edge_dir}
if [ ${REPLICAS} -eq 1 ]; then
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run1_protein/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run1_water/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run1_protein/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run1_water/pre
elif [ ${REPLICAS} -eq 2 ]; then
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run1_protein/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run2_protein/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run1_water/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run2_water/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run1_protein/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run2_protein/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run1_water/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run2_water/pre
elif [ ${REPLICAS} -eq 3 ]; then
        mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid
        mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run1_protein/pre
        mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run2_protein/pre
        mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run3_protein/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run1_water/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run2_water/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateA/run3_water/pre
        mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid
        mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run1_protein/pre
        mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run2_protein/pre
        mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run3_protein/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run1_water/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run2_water/pre
	mkdir -p ${TRANSITION_PATH}${edge_dir}/stateB/run3_water/pre
else
        echo "you should set the num of replicas."
	exit
fi	
done < ${edge_file}
else
echo "${edge_file} does not exist"
exit
fi


