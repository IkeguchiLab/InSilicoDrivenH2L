#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
    echo "\$PROJ_HOME/config.bash not found"
    exit
fi
source ${PROJ_HOME}config.bash

mkdir -p ${PROJ_HOME}workpath/job_prepare_transition
cp ${PROJ_HOME}/config.bash ${PROJ_HOME}workpath/job_prepare_transition/

edge_file=${TRANSITION_PATH}edgelist.txt
echo ${edge_file}

count=0
if [ -e "${edge_file}" ]; then
while IFS=',' read -r ligA ligB;do
edge_dir=edge_${ligA}_${ligB}
echo ${edge_dir}
echo ${count}

cat << END > ${PROJ_HOME}workpath/job_prepare_transition/jobscript_prepare_transition_${count}
#!/bin/bash
#SBATCH --job-name=${edge_dir}_prepare_transition
#SBATCH --get-user-env
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --partition=c4l

source config.bash
module load intelmpi

num_run=`find ${TRANSITION_PATH}${edge_dir}/stateA/ -type d -name "run*_protein" | wc -l`

source /share/local/miniforge3/etc/profile.d/conda.sh

conda activate AZtutorial

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/pairs1.dat" ]; then
pmx atomMapping -i1 ${LIGAND_PATH}${ligA}/mol_gmx.pdb -i2 ${LIGAND_PATH}${ligB}/mol_gmx.pdb -o1 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/pairs1.dat -o2 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/pairs2.dat -opdb1 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdb1.pdb -opdb2 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdb2.pdb -opdbm1 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdbm1.pdb -opdbm2 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdbm2.pdb -score ${TRANSITION_PATH}${edge_dir}/stateA/score.dat -log ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/mapping.log
echo "--d 0.05" > ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/mapping_condition
sleep 0.2
fi

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/pairs1.dat" ]; then
pmx atomMapping --d 0.075 -i1 ${LIGAND_PATH}${ligA}/mol_gmx.pdb -i2 ${LIGAND_PATH}${ligB}/mol_gmx.pdb -o1 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/pairs1.dat -o2 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/pairs2.dat -opdb1 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdb1.pdb -opdb2 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdb2.pdb -opdbm1 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdbm1.pdb -opdbm2 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdbm2.pdb -score ${TRANSITION_PATH}${edge_dir}/stateA/score.dat -log ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/mapping.log
echo "--d 0.075" > ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/mapping_condition
sleep 0.2
fi

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/pairs1.dat" ]; then
pmx atomMapping --d 0.1 -i1 ${LIGAND_PATH}${ligA}/mol_gmx.pdb -i2 ${LIGAND_PATH}${ligB}/mol_gmx.pdb -o1 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/pairs1.dat -o2 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/pairs2.dat -opdb1 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdb1.pdb -opdb2 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdb2.pdb -opdbm1 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdbm1.pdb -opdbm2 ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/out_pdbm2.pdb -score ${TRANSITION_PATH}${edge_dir}/stateA/score.dat -log ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/mapping.log
echo "--d 0.1" > ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/mapping_condition
sleep 0.2
fi

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/merged.itp" ]; then
pmx ligandHybrid -i1 ${LIGAND_PATH}${ligA}/mol_gmx.pdb -i2 ${LIGAND_PATH}${ligB}/mol_gmx.pdb -itp1 ${LIGAND_PATH}${ligA}/MOL.itp -itp2  ${LIGAND_PATH}${ligB}/MOL.itp -pairs ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/pairs1.dat -oA ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/mergedA.pdb -oB ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/mergedB.pdb -oitp ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/merged.itp -offitp ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/ffmerged.itp -log ${TRANSITION_PATH}${edge_dir}/stateA/ligandHybrid/hybrid.log
sleep 0.2
fi

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/pairs1.dat" ]; then
pmx atomMapping -i1 ${LIGAND_PATH}${ligB}/mol_gmx.pdb -i2 ${LIGAND_PATH}${ligA}/mol_gmx.pdb -o1 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/pairs1.dat -o2 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/pairs2.dat -opdb1 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdb1.pdb -opdb2 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdb2.pdb -opdbm1 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdbm1.pdb -opdbm2 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdbm2.pdb -score ${TRANSITION_PATH}${edge_dir}/stateB/score.dat -log ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/mapping.log
echo "--d 0.05" > ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/mapping_condition
sleep 0.2
fi

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/pairs1.dat" ]; then
pmx atomMapping --d 0.075 -i1 ${LIGAND_PATH}${ligB}/mol_gmx.pdb -i2 ${LIGAND_PATH}${ligA}/mol_gmx.pdb -o1 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/pairs1.dat -o2 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/pairs2.dat -opdb1 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdb1.pdb -opdb2 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdb2.pdb -opdbm1 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdbm1.pdb -opdbm2 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdbm2.pdb -score ${TRANSITION_PATH}${edge_dir}/stateB/score.dat -log ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/mapping.log
echo "--d 0.075" > ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/mapping_condition
sleep 0.2
fi

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/pairs1.dat" ]; then
pmx atomMapping --d 0.1 -i1 ${LIGAND_PATH}${ligB}/mol_gmx.pdb -i2 ${LIGAND_PATH}${ligA}/mol_gmx.pdb -o1 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/pairs1.dat -o2 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/pairs2.dat -opdb1 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdb1.pdb -opdb2 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdb2.pdb -opdbm1 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdbm1.pdb -opdbm2 ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/out_pdbm2.pdb -score ${TRANSITION_PATH}${edge_dir}/stateB/score.dat -log ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/mapping.log
echo "--d 0.1" > ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/mapping_condition
sleep 0.2
fi

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/merged.itp" ]; then
pmx ligandHybrid -i1 ${LIGAND_PATH}${ligB}/mol_gmx.pdb -i2 ${LIGAND_PATH}${ligA}/mol_gmx.pdb -itp1 ${LIGAND_PATH}${ligB}/MOL.itp -itp2  ${LIGAND_PATH}${ligA}/MOL.itp -pairs ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/pairs1.dat -oA ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/mergedA.pdb -oB ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/mergedB.pdb -oitp ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/merged.itp -offitp ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/ffmerged.itp -log ${TRANSITION_PATH}${edge_dir}/stateB/ligandHybrid/hybrid.log
sleep 0.2
fi

for ((i=1; i<=\$num_run; i+=1))
do
if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/frame80.gro" ]; then
\${GMX} trjconv -s ${EQ_PATH}${ligA}/run\$i/cal03_eq_npt_protein/tpr.tpr -f ${EQ_PATH}${ligA}/run\$i/cal03_eq_npt_protein/traj.trr -o ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/frame.gro -sep -ur compact -pbc mol -b 2250 << EOF >& ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/log_trjconv
System
EOF
sleep 0.2
mv ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/frame0.gro ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/frame80.gro
fi

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/frame80.gro" ]; then
${GMX} trjconv -s ${EQ_PATH}${ligA}/run\$i/cal03_eq_npt_water/tpr.tpr -f ${EQ_PATH}${ligA}/run\$i/cal03_eq_npt_water/traj.trr -o ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/frame.gro -sep -ur compact -pbc mol -b 2250 << EOF >& ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/log_trjconv
System
EOF
sleep 0.2
mv ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/frame0.gro ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/frame80.gro
fi

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/frame80.gro" ]; then
${GMX} trjconv -s ${EQ_PATH}${ligB}/run\$i/cal03_eq_npt_protein/tpr.tpr -f ${EQ_PATH}${ligB}/run\$i/cal03_eq_npt_protein/traj.trr -o ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/frame.gro -sep -ur compact -pbc mol -b 2250 << EOF >& ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/log_trjconv
System
EOF
sleep 0.2
mv ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/frame0.gro ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/frame80.gro
fi

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/frame80.gro" ]; then
${GMX} trjconv -s ${EQ_PATH}${ligB}/run\$i/cal03_eq_npt_water/tpr.tpr -f ${EQ_PATH}${ligB}/run\$i/cal03_eq_npt_water/traj.trr -o ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/frame.gro -sep -ur compact -pbc mol -b 2250 << EOF >& ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/log_trjconv
System
EOF
sleep 0.2
mv ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/frame0.gro ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/frame80.gro
fi
done
END
count=`expr $count + 1`
done < ${edge_file}
else
echo "${edge_file} does not exist"
exit
fi

num_edge=`cat ${edge_file} | wc -l`
cat << EOF > ${PROJ_HOME}workpath/job_prepare_transition/job_submit.sh
#!/bin/bash
for ((i=0; i<${num_edge}; i++)); do
        echo \$i
        sbatch jobscript_prepare_transition_\${i}
done
EOF


