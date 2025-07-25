# NES_AWS_slurm
This repository contains necessary input files, scripts, and environment settings to perform NES calculations using SLURM on AWS.

## Note
- The configuration and simulation scripts are tailored specifically for the **PDE9A** protein target.  
  If you plan to apply this workflow to a different target, you will need to modify the corresponding files, particularly those related to the protein and complex parameter for MD simulation setup.

- This workflow is designed for execution on **AWS EC2 spot instances** using **SLURM** job scheduling, tested with instance types such as `g4dn.xlarge` and `c4.xlarge`.

## Directory Overview

- **config.bash**  
  A configuration script to set directory paths and the number of NES runs (1-3).

- **ligands/**  
  Contains parameter files for the ligands used in NES calculations.

- **scripts/**  
  Includes scripts used to perform NES calculations.

- **mdppath/**  
  Contains MD simulation input files (`.mdp`) required for NES calculations.

- **protein_amber/**  
  Parameter files for the target protein (PDE9A) used in NES simulations.

- **conda_env/**  
  Lists of Python packages used in two Conda environments:  
  - One for generating ligand parameter files with AmberTools (AmberTools) 
  - One for running and analyzing NES calculations (AZtutorial)

## Generated Directory Structure

After running the setup and job submission scripts, the following directory structure will be created. This structure contains all input/output files for equilibrium and transition simulations, along with the corresponding SLURM job scripts.

```
ligands/
└── liglist.txt              # List of ligands used to generate parameter files
└── ligands/
    └── liglist.txt              # List of ligands used for equilibrium calculations

workpath/
├── eq/                          # Input/output files for equilibrium MD simulations  
├── transition/                  # Input/output files for transition simulations
└── edgelist.txt              # List of ligand pairs (separated by a comma) to be used in transition calculations
├── job_eq_script/              # SLURM job scripts for running equilibrium simulations
├── job_prepare_transition/     # SLURM job scripts to generate parameter files for transition simulations
├── job_prepare_transition_tpr/ # SLURM job scripts to generate input .tpr files for transition simulations
├── job_transition_script/      # SLURM job scripts for running transition simulations
└── summary/                    # Summary of results from transition simulations
```

### Description of Directories and Files

- **ligand/ligand/liglist.txt**  
  Defines the list of ligands to be used for equilibrium simulations.

- **workpath/transition/edgelist.txt**  
  Defines the list of ligands used for equilibrium calculations

- **workpath/eq/**  
  Stores input and output files for the equilibrium molecular dynamics (MD) simulations.

- **workpath/transition/**  
  Stores input and output files for the transition simulations.

- **workpath/job_eq_script/**  
  Contains SLURM job scripts to run equilibrium MD simulations.

- **workpath/job_prepare_transition/**  
  Contains SLURM job scripts to generate necessary parameter files for transition simulations.

- **workpath/job_prepare_transition_tpr/**  
  Contains SLURM job scripts to create `.tpr` input files for the transition simulations.

- **workpath/job_transition_script/**  
  Contains SLURM job scripts to execute the transition simulations.

- **workpath/summary/**  
  Stores summary results and analysis outputs from transition simulations.
  
## Usage

## Workflow Overview

The overall workflow consists of the following steps:

1. **Ligand Parameter Generation**  
   Generate topology and coordinate files for each ligand using AmberTools.

2. **Preparation of Complex and Water Systems**  
   Create topology files and MD input files for:
   - Protein–ligand complex
   - Ligand-only in water  
   These are used in equilibrium simulations.

3. **Equilibrium MD Simulations**  
   Run equilibrium molecular dynamics simulations for both the complex and ligand-in-water systems using SLURM.

4. **Non-equilibrium Transition Setup**  
   Generate parameter files and input files for transition (nonequilibrium) simulations in both complex and water environments.

5. **Transition MD Simulations**  
   Run non-equilibrium transition MD simulations using SLURM.

6. **Summary and Analysis**  
   Summarize the results of the transition simulations to evaluate relative binding free energies.

## Conda Environments

This workflow requires two separate Conda environments:

### AmberTools Environment
Used for generating ligand parameter files from `.pdb` structures using `antechamber` and related AmberTools utilities.

Example environment file: 
conda_env/ambertools_env.yaml

### NES Calculation and Analysis Environment (AZtutorial)
Used to perform NES calculations and analyze the results using pmx.
For further information on how to use pmx, please refer to the official repository: 
🔗 https://github.com/deGrootLab/pmx/tree/master

Example environment file: 
conda_env/AZtutorial_env.yaml

## Step-by-Step Instructions

### 0: Configuration File Setup
Before starting the workflow, you must configure the paths and settings used throughout the simulations.

Edit the `config.bash` file to match your environment. Below is an example:

```bash
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

REPLICAS=3  # the number of NES runs
```
### 1: Ligand Parameter Generation
### 1-1. Set Environment Variables
``` bash
source config.bash
```

### 1-2. Prepare Ligand PDB Files

Place all ligand `.pdb` files into the `ligands/` directory.

### 1-3. Generate `liglist.txt`

Run the provided script to create `liglist.txt`, which lists all ligands to be used in generating ligand parameter files.

```bash
cd ligands
chmod +x ./create_liglist.sh
./create_liglist.sh
```

### 1-4. Generate Ligand Parameter Files

Submit the job to SLURM to generate ligand parameter files.

```bash
chmod +x ./prepare_lig.sh
sbatch prepare_lig.job
```

### 1-5. Check and Clean Ligand Parameter Files

Before proceeding to Step 2, check whether all ligand parameter files were successfully created.

```bash
cd ligands
chmod +x ./check_empty_mol_itp.sh
./check_empty_mol_itp.sh # This script identifies ligands for which .itp files were not properly generated
python itp_MOL_correction.py # edit .itp files to remove [ mol ] section
```

Next, copy the ligand list and manually update it.

```bash
cp ../liglist.txt ./
vi liglist.txt  # Remove ligands that failed parameter generation
```

### 2: Preparation of Complex and Water Systems
### 2-1. Create directories

```bash
cd ${PROJ_HOME}
chmod +x ./scripts/step00_directory.sh
./scripts/step00_directory.sh
```

### 2-2. Generate Parameter Files and Input Files for MD simulation (Equilibrium)

```bash
chmod +x ./scripts/step01_topology.sh
./scripts/step01_topology.sh
chmod +x ./scripts/step02_prepare_simulation.sh
./scripts/step02_prepare_simulation.sh
chmod +x ./scripts/step03_prepare_em.sh
./scripts/step03_prepare_em.sh
```

### 3. Equilibrium MD Simulations

```bash
chmod +x ./scripts/step04_jobscripts_eq.sh
./scripts/step04_jobscripts_eq.sh
cd workpath/job_eq_script
chmod +x ./job_submit.sh
./job_submit.sh
```

### 4. Non-equilibrium Transition Setup
### 4-1. Create lists for the transition simulations
Before generating input files for the transition simulations, define the ligand transformation network by creating a list of ligand pairs.

Create or edit the following file:

```bash
workpath/transition/edgelist.txt
```

Each line should specify a ligand pair (comma-separated) to be used in transition (nonequilibrium) simulations. For example:

```bash
ligA,ligB
ligB,ligC
ligA,ligC
```

### 4-2. Create directories

```bash
cd {PROJ_HOME}
chmod +x scripts/step05_transition_directory.sh
scripts/step05_transition_directory.sh
```

### 4-3. Generate Hybrid Topology Files

```bash
chmod +x scripts/step06_prepare_transition_per_edge.sh
scripts/step06_prepare_transition_per_edge.sh
cd workpath/job_prepare_transition
chmod +x ./job_submit.sh
./job_submit.sh
```

### 4-4. Generate Input Files for the Transition Simulations

```bash
cd ${PROJ_HOME}
chmod +x scripts/step07_prepare_transition_jobscript_per_edge.sh
scripts/step07_prepare_transition_jobscript_per_edge.sh
cd workpath/job_prepare_transition_tpr
chmod +x ./job_submit.sh
./job_submit.sh
```

### 5. Transition MD Simulations

```bash
cd ${PROJ_HOME}
chmod +x scripts/step08_jobscripts_transition.sh
scripts/step08_jobscripts_transition.sh
cd workpath/job_transition_script
chmod +x ./job_submit.sh
./job_submit.sh
```

### 6. Summary and Analysis

```bash
cd ${PROJ_HOME}
chmod +x scripts/step09_pmx_analyse.sh
scripts/step09_pmx_analyse.sh
chmod +x scripts/step10_summary_result.sh
scripts/step10_summary_result
chmod +x scripts/step11_summary_analyse.sh
scripts/step11_summary_analyse
```

This script will generate the final output file.

```bash
workpath/summary/results.csv
```
