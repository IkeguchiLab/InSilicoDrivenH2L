#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
    echo "\$PROJ_HOME/config.bash not found"
    exit
fi
source ${PROJ_HOME}config.bash

protein_em_tpr=`find ${PROJ_HOME}workpath/eq/lig_*/run*/cal01_em_protein/ -type f -name "tpr.tpr"`
protein_em_tpr_array=(${protein_em_tpr})
len_protein_em_tpr_array=${#protein_em_tpr_array[@]}

water_em_tpr=`find ${PROJ_HOME}workpath/eq/lig_*/run*/cal01_em_water/ -type f -name "tpr.tpr"`
water_em_tpr_array=(${water_em_tpr})
len_water_em_tpr_array=${#water_em_tpr_array[@]}

mkdir -p ${PROJ_HOME}workpath/job_eq_script
cp ${PROJ_HOME}/config.bash ${PROJ_HOME}workpath/job_eq_script/

for ((i=0; i<len_protein_em_tpr_array; i++));do
echo ${protein_em_tpr_array[i]}
protein_em_dir=$(dirname "${protein_em_tpr_array[i]}")
echo $protein_em_dir

cat << EOF > ${PROJ_HOME}workpath/job_eq_script/jobscript_eq_protein_${i}
#!/bin/bash
#SBATCH --job-name=${protein_em_dir}
#SBATCH --get-user-env
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --partition=g4dn

source ./config.bash
module load intelmpi

export GMXRUN="$GMX mdrun"
log_file="md.log"
tpr_file="tpr.tpr"
cd ${protein_em_dir}
if [ -f "\${log_file}" ];then
    if ! tail -n 2 "\${log_file}" | grep -q "Finished";then
        \$GMXRUN -s tpr.tpr
    else
        echo "log file exists."
    fi
else
    \$GMXRUN -s tpr.tpr
fi
rm \#*

cd ../

mkdir cal02_eq_nvt_protein
cd cal02_eq_nvt_protein
if [ -f "\${log_file}" ];then
    if ! tail -n 2 "\${log_file}" | grep -q "Finished";then
        \$GMXRUN -s tpr.tpr
    else
        echo "log file exists."
    fi
else
    ${GMX} grompp -f ${MDP_PATH}eq_nvt_l0.mdp -c ../cal01_em_protein/confout.gro -r ../cal01_em_protein/confout.gro -p ../../protein/topol.top -o tpr.tpr -maxwarn 4 -po mdout.mdp
    sleep 0.2
    \$GMXRUN -s tpr.tpr
fi
rm \#*

cd ../

mkdir cal03_eq_npt_protein
cd cal03_eq_npt_protein
if [ -f "\${log_file}" ];then
    if ! tail -n 2 "\${log_file}" | grep -q "Finished";then
        shopt -s nullglob
        check_files=(*.cpt)

        if [[ \${#check_files[@]} -eq 0 ]]; then
            \$GMXRUN -s tpr.tpr -cpi \$SLURM_JOBID.cpt -cpo \$SLURM_JOBID.cpt
        else
            for check_file in *.cpt; do
                if [[ \${check_file} != *_prev.cpt ]]; then
                    \$GMXRUN -s tpr.tpr -cpi \${check_file} -cpo \${check_file}
                fi
            done
        fi
    else
        echo "log file exists."
    fi
else
    ${GMX} grompp -f ${MDP_PATH}eq_npt_l0.mdp -c ../cal02_eq_nvt_protein/confout.gro -p ../../protein/topol.top -o tpr.tpr -maxwarn 4 -po mdout.mdp
    sleep 0.2
    \$GMXRUN -s tpr.tpr -cpi \$SLURM_JOBID.cpt -cpo \$SLURM_JOBID.cpt
fi
rm \#

cd ../
EOF
done

for ((i=0; i<len_water_em_tpr_array; i++));do
echo ${water_em_tpr_array[i]}
water_em_dir=$(dirname "${water_em_tpr_array[i]}")
echo ${water_em_dir}

cat << EOF > ${PROJ_HOME}workpath/job_eq_script/jobscript_eq_water_${i}
#!/bin/bash
#SBATCH --job-name=${water_em_dir}
#SBATCH --get-user-env
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --partition=g4dn

source ./config.bash
module load intelmpi

export GMXRUN="$GMX mdrun"
cd ${water_em_dir}
\$GMXRUN -s tpr.tpr

cd ../

mkdir cal02_eq_nvt_water
cd cal02_eq_nvt_water
${GMX} grompp -f ${MDP_PATH}eq_nvt_l0.mdp -c ../cal01_em_water/confout.gro -r ../cal01_em_water/confout.gro -p ../../water/topol.top -o tpr.tpr -maxwarn 4 -po mdout.mdp
sleep 0.2
\$GMXRUN -s tpr.tpr

cd ../

mkdir cal03_eq_npt_water
cd cal03_eq_npt_water
${GMX} grompp -f ${MDP_PATH}eq_npt_l0.mdp -c ../cal02_eq_nvt_water/confout.gro -p ../../water/topol.top -o tpr.tpr -maxwarn 4 -po mdout.mdp
sleep 0.2
\$GMXRUN -s tpr.tpr

cd ../
EOF
done

cat << EOF > ${PROJ_HOME}workpath/job_eq_script/job_submit.sh
#!/bin/bash
for ((i=0; i<${len_protein_em_tpr_array}; i++)); do
        echo \$i
        sbatch jobscript_eq_protein_\${i}
        sbatch jobscript_eq_water_\${i}
done
EOF
