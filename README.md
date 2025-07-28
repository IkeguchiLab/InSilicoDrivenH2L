# **In Silico-Driven H2L: A Case Study on PDE9A Inhibitors**
## Overview

This repository contains:

- Scripts for Non-equilibrium Switching (NES) calculations  
- A list of selected R-groups  
- NES results obtained at steps 1-3  
- Scripts used to generate compounds at steps 1-3

## Directory Structure
```
InSilicoDrivenH2L/
├── NES_AWS_slurm/ # Scripts for running NES on AWS cloud using Slurm
├── R_group/ # Selected 30,477 R-groups
├── RBFE_validation/ # Performance evaluation results of NES, FEP+, and OpenFE
├── Step1/ # Compound generation on Cores 1-4 and NES results of representative compounds
├── Step2/ # H pocket optimization and screening results
├── Step3/ # M pocket optimization and screening results
├── conda_env/ # Conda virtual environment used for the study
├── Appendix/ # tautomer, UGT
└── README.md
```
