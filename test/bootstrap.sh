#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=2
cd "/home/jsicherman/Thesis Work"

Rscript "test/bootstrap.R"
