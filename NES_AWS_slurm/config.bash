#!/bin/bash

PROJ_HOME="/path/to/your/directory/"
WORK_PATH="${PROJ_HOME}workpath/"
MDP_PATH="${PROJ_HOME}mdppath/"

GMX="/share/local/gromacs/gromacs-2024.1-impi/bin/gmx"

PROTEIN_PATH="${PROJ_HOME}protein_amber/"
LIGAND_PATH="${PROJ_HOME}ligands/ligands/"
EQ_PATH="${PROJ_HOME}eq/"
TRANSITION_PATH="${PROJ_HOME}workpath/transition/"
SCRIPTS_PATH="${PROJ_HOME}scripts/"
REPLICAS=3 #the number of runs
