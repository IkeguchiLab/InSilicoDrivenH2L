#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
    echo "\$PROJ_HOME/config.bash not found"
    exit
fi
source ${PROJ_HOME}config.bash

mkdir -p ${PROJ_HOME}workpath/job_prepare_transition_tpr
cp ${PROJ_HOME}/config.bash ${PROJ_HOME}workpath/job_prepare_transition_tpr/

edge_file=${TRANSITION_PATH}edgelist.txt
echo ${edge_file}

count=0
if [ -e "${edge_file}" ]; then
while IFS=',' read -r ligA ligB;do
edge_dir=edge_${ligA}_${ligB}
echo ${edge_dir}
echo ${count}

cat << END > ${PROJ_HOME}workpath/job_prepare_transition_tpr/jobscript_prepare_transition_tpr_${count}
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

for ((i=1; i<=\$num_run; i+=1))
do

###
# transition protein A to B 
###

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/ti80.tpr" ]; then
cd ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/
cp ${EQ_PATH}${ligA}/protein/topol.top ./
ffmerged_itp_path="${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/ffmerged.itp"
merged_itp_path="${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/merged.itp"
sed -i "2s|.*|#include \"\${ffmerged_itp_path}\"|" "topol.top"
sed -i "3s|.*|#include \"\${merged_itp_path}\"|" "topol.top"

for i2 in {1..80};do
    gmx editconf -f frame\${i2}.gro -o frame\${i2}.pdb
    cat frame\${i2}.pdb | grep MOL > frame\${i2}_${ligA}.pdb
    gmx editconf -f frame\${i2}_${ligA}.pdb -o frame\${i2}_${ligA}_edit.pdb
    sleep 0.02
    pmx ligandHybrid -i1 frame\${i2}_${ligA}_edit.pdb -i2 ${LIGAND_PATH}${ligB}/mol_gmx.pdb -itp1 ${LIGAND_PATH}${ligA}/MOL.itp -itp2 ${LIGAND_PATH}${ligB}/MOL.itp -pairs ../../ligandHybrid/pairs1.dat -oA mergedA_frame\${i2}.pdb -oB mergedB_frame\${i2}.pdb -oitp merged.itp -offitp ffmerged_pmx.itp
    
    cat ${LIGAND_PATH}${ligA}/ffMOL.itp ${LIGAND_PATH}${ligB}/ffMOL.itp | awk '!a[\$0]++{print}' > ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/ffMOL_${ligA}_${ligB}.itp
    sed -i '1d' "ffmerged_pmx.itp"
    cat ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/ffMOL_${ligA}_${ligB}.itp ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/ffmerged_pmx.itp > ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_protein/pre/ffmerged.itp

    cat mergedA_frame\${i2}.pdb | grep MOL > mergedA_frame\${i2}_edit.pdb
    python ${SCRIPTS_PATH}edit_MOL.py frame\${i2}.pdb frame\${i2}_${ligA}_merged.pdb mergedA_frame\${i2}_edit.pdb
    gmx editconf -f frame\${i2}_${ligA}_merged.pdb -o frame\${i2}_${ligA}_merged.gro
    gmx grompp -f ${MDP_PATH}transitions_from_npt_l0_05fs.mdp -c frame\${i2}_${ligA}_merged.gro -r frame\${i2}_${ligA}_merged.gro -p topol.top -o ../ti\${i2}.tpr -po ../mdout_\${i2}.mdp -maxwarn 7
done
rm \#*
rm ../\#*
rm *.pdb
rm ../*.pdb
rm *_merged.gro
fi

###
# transition protein B to A
###

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/ti80.tpr" ]; then
cd ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/
cp ${EQ_PATH}${ligB}/protein/topol.top ./
ffmerged_itp_path="${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/ffmerged.itp"
merged_itp_path="${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/merged.itp"
sed -i "2s|.*|#include \"\${ffmerged_itp_path}\"|" "topol.top"
sed -i "3s|.*|#include \"\${merged_itp_path}\"|" "topol.top"

for i2 in {1..80};do
    gmx editconf -f frame\${i2}.gro -o frame\${i2}.pdb
    cat frame\${i2}.pdb | grep MOL > frame\${i2}_${ligB}.pdb
    gmx editconf -f frame\${i2}_${ligB}.pdb -o frame\${i2}_${ligB}_edit.pdb
    sleep 0.02
    pmx ligandHybrid -i1 frame\${i2}_${ligB}_edit.pdb -i2 ${LIGAND_PATH}${ligA}/mol_gmx.pdb -itp1 ${LIGAND_PATH}${ligB}/MOL.itp -itp2 ${LIGAND_PATH}${ligA}/MOL.itp -pairs ../../ligandHybrid/pairs1.dat -oA mergedA_frame\${i2}.pdb -oB mergedB_frame\${i2}.pdb -oitp merged.itp -offitp ffmerged_pmx.itp

    cat ${LIGAND_PATH}${ligB}/ffMOL.itp ${LIGAND_PATH}${ligA}/ffMOL.itp | awk '!a[\$0]++{print}' > ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/ffMOL_${ligB}_${ligA}.itp
    sed -i '1d' "ffmerged_pmx.itp"
    cat ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/ffMOL_${ligB}_${ligA}.itp ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/ffmerged_pmx.itp > ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_protein/pre/ffmerged.itp

    cat mergedA_frame\${i2}.pdb | grep MOL > mergedA_frame\${i2}_edit.pdb
    python ${SCRIPTS_PATH}edit_MOL.py frame\${i2}.pdb frame\${i2}_${ligB}_merged.pdb mergedA_frame\${i2}_edit.pdb
    gmx editconf -f frame\${i2}_${ligB}_merged.pdb -o frame\${i2}_${ligB}_merged.gro
    gmx grompp -f ${MDP_PATH}transitions_from_npt_l0_05fs.mdp -c frame\${i2}_${ligB}_merged.gro -r frame\${i2}_${ligB}_merged.gro -p topol.top -o ../ti\${i2}.tpr -po ../mdout_\${i2}.mdp -maxwarn 7
done
rm \#*
rm ../\#*
rm *.pdb
rm ../*.pdb
rm *_merged.gro
fi

###
# transition water A to B
###

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/ti80.tpr" ]; then
cd ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/
cp ${EQ_PATH}${ligA}/water/topol.top ./
ffmerged_itp_path="${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/ffmerged.itp"
merged_itp_path="${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/merged.itp"
sed -i "2s|.*|#include \"\${ffmerged_itp_path}\"|" "topol.top"
sed -i "3s|.*|#include \"\${merged_itp_path}\"|" "topol.top"

for i2 in {1..80};do
    gmx editconf -f frame\${i2}.gro -o frame\${i2}.pdb
    cat frame\${i2}.pdb | grep MOL > frame\${i2}_${ligA}.pdb
    gmx editconf -f frame\${i2}_${ligA}.pdb -o frame\${i2}_${ligA}_edit.pdb
    sleep 0.02
    pmx ligandHybrid -i1 frame\${i2}_${ligA}_edit.pdb -i2 ${LIGAND_PATH}${ligB}/mol_gmx.pdb -itp1 ${LIGAND_PATH}${ligA}/MOL.itp -itp2 ${LIGAND_PATH}${ligB}/MOL.itp -pairs ../../ligandHybrid/pairs1.dat -oA mergedA_frame\${i2}.pdb -oB mergedB_frame\${i2}.pdb -oitp merged.itp -offitp ffmerged_pmx.itp

    cat ${LIGAND_PATH}${ligA}/ffMOL.itp ${LIGAND_PATH}${ligB}/ffMOL.itp | awk '!a[\$0]++{print}' > ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/ffMOL_${ligA}_${ligB}.itp
    sed -i '1d' "ffmerged_pmx.itp"
    cat ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/ffMOL_${ligA}_${ligB}.itp ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/ffmerged_pmx.itp > ${TRANSITION_PATH}${edge_dir}/stateA/run\${i}_water/pre/ffmerged.itp

    cat mergedA_frame\${i2}.pdb | grep MOL > mergedA_frame\${i2}_edit.pdb
    python ${SCRIPTS_PATH}edit_MOL.py frame\${i2}.pdb frame\${i2}_${ligA}_merged.pdb mergedA_frame\${i2}_edit.pdb
    gmx editconf -f frame\${i2}_${ligA}_merged.pdb -o frame\${i2}_${ligA}_merged.gro
    gmx grompp -f ${MDP_PATH}transitions_from_npt_l0_05fs.mdp -c frame\${i2}_${ligA}_merged.gro -r frame\${i2}_${ligA}_merged.gro -p topol.top -o ../ti\${i2}.tpr -po ../mdout_\${i2}.mdp -maxwarn 7
done
rm \#*
rm ../\#*
rm *.pdb
rm ../*.pdb
rm *_merged.gro
fi

###
# transition water B to A
###

if [ ! -e "${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/ti80.tpr" ]; then
cd ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/
cp ${EQ_PATH}${ligB}/water/topol.top ./
ffmerged_itp_path="${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/ffmerged.itp"
merged_itp_path="${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/merged.itp"
sed -i "2s|.*|#include \"\${ffmerged_itp_path}\"|" "topol.top"
sed -i "3s|.*|#include \"\${merged_itp_path}\"|" "topol.top"

for i2 in {1..80};do
    gmx editconf -f frame\${i2}.gro -o frame\${i2}.pdb
    cat frame\${i2}.pdb | grep MOL > frame\${i2}_${ligB}.pdb
    gmx editconf -f frame\${i2}_${ligB}.pdb -o frame\${i2}_${ligB}_edit.pdb
    sleep 0.02
    pmx ligandHybrid -i1 frame\${i2}_${ligB}_edit.pdb -i2 ${LIGAND_PATH}${ligA}/mol_gmx.pdb -itp1 ${LIGAND_PATH}${ligB}/MOL.itp -itp2 ${LIGAND_PATH}${ligA}/MOL.itp -pairs ../../ligandHybrid/pairs1.dat -oA mergedA_frame\${i2}.pdb -oB mergedB_frame\${i2}.pdb -oitp merged.itp -offitp ffmerged_pmx.itp
    cat mergedA_frame\${i2}.pdb | grep MOL > mergedA_frame\${i2}_edit.pdb

    cat ${LIGAND_PATH}${ligB}/ffMOL.itp ${LIGAND_PATH}${ligA}/ffMOL.itp | awk '!a[\$0]++{print}' > ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/ffMOL_${ligB}_${ligA}.itp
    sed -i '1d' "ffmerged_pmx.itp"
    cat ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/ffMOL_${ligB}_${ligA}.itp ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/ffmerged_pmx.itp > ${TRANSITION_PATH}${edge_dir}/stateB/run\${i}_water/pre/ffmerged.itp

    python ${SCRIPTS_PATH}edit_MOL.py frame\${i2}.pdb frame\${i2}_${ligB}_merged.pdb mergedA_frame\${i2}_edit.pdb
    gmx editconf -f frame\${i2}_${ligB}_merged.pdb -o frame\${i2}_${ligB}_merged.gro
    gmx grompp -f ${MDP_PATH}transitions_from_npt_l0_05fs.mdp -c frame\${i2}_${ligB}_merged.gro -r frame\${i2}_${ligB}_merged.gro -p topol.top -o ../ti\${i2}.tpr -po ../mdout_\${i2}.mdp -maxwarn 7
done
rm \#*
rm ../\#*
rm *.pdb
rm ../*.pdb
rm *_merged.gro
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
cat << EOF > ${PROJ_HOME}workpath/job_prepare_transition_tpr/job_submit.sh
#!/bin/bash
for ((i=0; i<${num_edge}; i++)); do
        echo \$i
        sbatch jobscript_prepare_transition_tpr_\${i}
done
EOF









