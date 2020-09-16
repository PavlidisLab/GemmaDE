#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=30
cd "/home/jsicherman/Thesis Work"

Rscript "test/artificial.R"
