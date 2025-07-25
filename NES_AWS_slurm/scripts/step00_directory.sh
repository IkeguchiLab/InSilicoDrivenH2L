#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
    echo "\$PROJ_HOME/config.bash not found"
    exit
fi
source ${PROJ_HOME}config.bash

if [ -f "${LIGAND_PATH}liglist.txt" ]; then
    while read lig
    do
    if [ ! -d ${EQ_PATH}${lig} ]; then
        mkdir -p ${EQ_PATH}${lig}
    else
        echo "${EQ_PATH}${lig} is already exist"
    fi
    done < "${LIGAND_PATH}liglist.txt"
else
    echo "${LIGAND_PATH}liglist.txt is not exist"
    exit
fi

if [ ${REPLICAS} -eq 1 ]; then
    while read lig
    do
        mkdir -p ${EQ_PATH}${lig}/run1/
    done < "${LIGAND_PATH}liglist.txt"
elif [ ${REPLICAS} -eq 2 ]; then
    while read lig
    do
        mkdir -p ${EQ_PATH}${lig}/run1/
        mkdir -p ${EQ_PATH}${lig}/run2/
    done < "${LIGAND_PATH}liglist.txt"
elif [ ${REPLICAS} -eq 3 ]; then
    while read lig
    do
        mkdir -p ${EQ_PATH}${lig}/run1/
        mkdir -p ${EQ_PATH}${lig}/run2/
        mkdir -p ${EQ_PATH}${lig}/run3/
    done < "${LIGAND_PATH}liglist.txt"
else
    echo "you should set the num of replicas."
fi
	
