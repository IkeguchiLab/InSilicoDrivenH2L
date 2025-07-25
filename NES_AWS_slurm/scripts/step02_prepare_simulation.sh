#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
    echo "\$PROJ_HOME/config.bash not found"
    exit
fi
source ${PROJ_HOME}config.bash

while read lig
do

###
# prepare simulation of complex
###

${GMX} editconf -f ${EQ_PATH}${lig}/protein/init.pdb -o ${EQ_PATH}${lig}/protein/box.pdb -bt dodecahedron -d 1.5 >& ${EQ_PATH}${lig}/protein/log_editconf
sleep 0.2

${GMX} solvate -cp ${EQ_PATH}${lig}/protein/box.pdb -cs spc216.gro -p ${EQ_PATH}${lig}/protein/topol.top -o ${EQ_PATH}${lig}/protein/water.pdb >& ${EQ_PATH}${lig}/protein/log_solvate
sleep 0.2

${GMX} grompp -f ${MDP_PATH}em_l0.mdp -c ${EQ_PATH}${lig}/protein/water.pdb -r ${EQ_PATH}${lig}/protein/water.pdb -p ${EQ_PATH}${lig}/protein/topol.top -o ${EQ_PATH}${lig}/protein/tpr.tpr -maxwarn 4 -po ${EQ_PATH}${lig}/protein/mdput.mdp >& ${EQ_PATH}${lig}/protein/log_grompp_ion
sleep 0.2

${GMX} genion -s ${EQ_PATH}${lig}/protein/tpr.tpr -p ${EQ_PATH}${lig}/protein/topol.top -o ${EQ_PATH}${lig}/protein/ions.pdb -conc 0.15 -neutral -pname NaJ -nname ClJ << EOF >& ${EQ_PATH}${lig}/protein/log_genion
SOL
EOF

rm ${EQ_PATH}${lig}/protein/\#*

###
# prepare simulation of ligand
###

${GMX} editconf -f ${EQ_PATH}${lig}/water/init.pdb -o ${EQ_PATH}${lig}/water/box.pdb -bt dodecahedron -d 1.5 >& ${EQ_PATH}${lig}/water/log_editconf
sleep 0.2

${GMX} solvate -cp ${EQ_PATH}${lig}/water/box.pdb -cs spc216.gro -p ${EQ_PATH}${lig}/water/topol.top -o ${EQ_PATH}${lig}/water/water.pdb >& ${EQ_PATH}${lig}/water/log_solvate >&${EQ_PATH}${lig}/water/log_solvate
sleep 0.2

${GMX} grompp -f ${MDP_PATH}em_l0.mdp -c ${EQ_PATH}${lig}/water/water.pdb -r ${EQ_PATH}${lig}/water/water.pdb -p ${EQ_PATH}${lig}/water/topol.top -o ${EQ_PATH}${lig}/water/tpr.tpr -maxwarn 4 -po ${EQ_PATH}${lig}/water/mdput.mdp >& ${EQ_PATH}${lig}/water/log_grompp_ion
sleep 0.2

${GMX} genion -s ${EQ_PATH}${lig}/water/tpr.tpr -p ${EQ_PATH}${lig}/water/topol.top -o ${EQ_PATH}${lig}/water/ions.pdb -conc 0.15 -neutral -pname NaJ -nname ClJ << EOF >& ${EQ_PATH}${lig}/water/log_genion
SOL
EOF

rm ${EQ_PATH}${lig}/water/\#*

done < "${LIGAND_PATH}liglist.txt"


