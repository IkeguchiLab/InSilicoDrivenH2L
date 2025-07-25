#!/bin/bash

if [ ! -f ${PROJ_HOME}config.bash ]; then
echo "\$PROJ_HOME/config.bash not found"
exit
fi
source ${PROJ_HOME}config.bash

mkdir -p ${PROJ_HOME}workpath/job_transition_script

cp ${PROJ_HOME}/config.bash ${PROJ_HOME}workpath/job_transition_script/

transition_dir=`find ${TRANSITION_PATH}edge*/state*/run*_* -type d -name "run*_*"`
transition_dir_array=(${transition_dir})
len_transition_dir_array=${#transition_dir_array[@]}

for ((i=0; i<len_transition_dir_array; i++));do
#echo ${transition_dir_array[i]}

dir1=$(basename "${transition_dir_array[i]}")
dir2=$(basename $(dirname "${transition_dir_array[i]}"))
dir3=$(basename $(dirname $(dirname "${transition_dir_array[i]}")))

formatted_var="${dir3}_${dir2}_${dir1}"

cat << EOF > ${PROJ_HOME}workpath/job_transition_script/jobscript_transition_${i}
#!/bin/bash
#SBATCH --job-name=${transition_dir_array[i]}
#SBATCH --get-user-env
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH --partition=g4dn

#SBATCH --gres=gpu:1

source config.bash
module load intelmpi

export GMXRUN="\$GMX mdrun"
cd ${transition_dir_array[i]}
for i in {1..80};do
dhdl_file="dhdl\${i}.xvg"
if [ ! -e \${dhdl_file} ]; then
\$GMXRUN -s ti\${i}.tpr -dhdl dhdl\${i} -g md\${i}.log
else
last_line=\$(tail -n 1 "\${dhdl_file}")
first_two_charts="\${last_line:0:2}"
if [ "\${first_two_charts}" = "50" ];then
echo "\${dhdl_file} is already exist."
else
\$GMXRUN -s ti\${i}.tpr -dhdl dhdl\${i} -g md\${i}.log
fi
fi
done
rm \#*
EOF
done

cat << EOF > ${PROJ_HOME}workpath/job_transition_script/job_submit.sh
#!/bin/bash
for ((i=0; i<${len_transition_dir_array}; i++)); do
        echo \$i
        sbatch jobscript_transition_\${i}
done
EOF
