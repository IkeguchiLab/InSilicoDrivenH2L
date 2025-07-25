#!/bin/bash

dir_path=$(pwd)

while read lig
do
    mkdir -p ligands/${lig}/pre
    mv ${lig}.pdb ligands/${lig}/pre
    cd ligands/${lig}/pre
    cat ${lig}.pdb | grep "HETATM" > ${lig}_grep.pdb
    charge=$(python ${dir_path}/../scripts/mol_charge.py ${dir_path}/ligands/${lig}/pre/${lig}.pdb)
    antechamber -fi pdb -fo prepi -i ${lig}_grep.pdb -o ${lig}.prep -c bcc -at gaff2 -nc ${charge} > antechamber_log
    parmchk2 -i ${lig}.prep -o ${lig}.frcmod -f prepi -s gaff2
    cat << EOF > tleap.in
    source leaprc.gaff2
    loadamberprep ${lig}.prep
    loadamberparams ${lig}.frcmod
    mol = loadpdb NEWPDB.PDB
    saveamberparm mol ${lig}.prmtop ${lig}.inpcrd
    quit
EOF
    tleap -f tleap.in
    python ${dir_path}/../scripts/parmed_2gro.py ${lig}
    sed -e "s/UNK/MOL/g" gromacs.pdb > mol_gmx.pdb
    mv mol_gmx.pdb ../
    sed -e "s/UNK/MOL/g" gromacs.top > gromacs_sed.top
    python ${dir_path}/../scripts/itp_from_gmx.py
    mv MOL.itp ffMOL.itp ../
    cd ../../../
done < liglist.txt
