#!/bin/bash
# This file can be used to pass in a variable path to run all generate scripts
# Rscript generate/fbm.R
Rscript generate/compile.R >generate/compile_log 2>generate/compile_log_err
Rscript generate/bootstrap.R >generate/bootstrap_log 2>generate/bootstrap_log_err
Rscript generate/package_data.R >generate/package_data_log 2>generate/package_data_log_err
