#!/bin/bash
# This file can be used to pass in a variable path to run all generate scripts
Rscript generate/fbm.R
Rscript generate/bootstrap.R
Rscript generate/package_data.R
