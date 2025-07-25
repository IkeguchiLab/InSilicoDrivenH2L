#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
    echo "\$PROJ_HOME/config.bash not found"
    exit
fi
source ${PROJ_HOME}config.bash

while read lig
do
mkdir -p ${EQ_PATH}${lig}/protein/
mkdir -p ${EQ_PATH}${lig}/water/

###
# topology of complex
###

cat << EOF > ${EQ_PATH}${lig}/protein/topol.top
#include "amber99sb-star-ildn-mut.ff/forcefield.itp"
#include "${LIGAND_PATH}${lig}/ffMOL.itp"
#include "${LIGAND_PATH}${lig}/MOL.itp"
#include "${PROTEIN_PATH}topol_Protein_chain_A.itp"
#include "${PROTEIN_PATH}topol_Ion_chain_A2.itp"
#include "${PROTEIN_PATH}topol_Other_chain_A3.itp"
#include "amber99sb-star-ildn-mut.ff/tip3p.itp"
#include "amber99sb-star-ildn-mut.ff/ions.itp"

[ system ]
protein and ligand in water

[ molecules ]
Protein_chain_A 1
Ion_chain_A2 1
Other_chain_A3 1
MOL 1
SOL 4
EOF

###
# topology of ligand
###

cat << EOF > ${EQ_PATH}${lig}/water/topol.top
#include "amber99sb-star-ildn-mut.ff/forcefield.itp"
#include "${LIGAND_PATH}${lig}/ffMOL.itp"
#include "${LIGAND_PATH}${lig}/MOL.itp"
#include "amber99sb-star-ildn-mut.ff/tip3p.itp"
#include "amber99sb-star-ildn-mut.ff/ions.itp"

[ system ]
ligand in water

[ molecules ]
MOL 1
EOF

###
# cordinate of complex
###

cat ${PROTEIN_PATH}protein.pdb ${LIGAND_PATH}${lig}/mol_gmx.pdb | grep -E 'ATOM|HETATM' > ${EQ_PATH}${lig}/protein/init_conc.pdb
cat ${EQ_PATH}${lig}/protein/init_conc.pdb | grep -v "HOH" > ${EQ_PATH}${lig}/protein/init_notHOH.pdb
cat ${EQ_PATH}${lig}/protein/init_conc.pdb | grep "HOH" > ${EQ_PATH}${lig}/protein/HOH.pdb
cat ${EQ_PATH}${lig}/protein/init_notHOH.pdb ${EQ_PATH}${lig}/protein/HOH.pdb > ${EQ_PATH}${lig}/protein/init_HOH.pdb
sed "s/HOH/SOL/g" ${EQ_PATH}${lig}/protein/init_HOH.pdb > ${EQ_PATH}${lig}/protein/init.pdb

###
# cordinate of ligand
###

cp ${LIGAND_PATH}${lig}/mol_gmx.pdb ${EQ_PATH}${lig}/water/init.pdb
done < "${LIGAND_PATH}liglist.txt"
