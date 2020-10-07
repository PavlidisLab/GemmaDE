#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=12
cd "/home/jsicherman/Thesis Work"

Rscript "test/artificial2.R"
