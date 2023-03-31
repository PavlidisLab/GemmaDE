#!/bin/bash
# This file can be used to pass in a variable path to run all generate scripts
# Rscript generate/fbm.R # old method to generate DATA.HOLDER. no longer used
Rscript generate/compile.R >generate/compile_log 2>generate/compile_log_err # complile the data from Gemma using the API
Rscript generate/term_fixes.R >generate/term_fixes 2>generate/term_fixes_err # fix ontology terms that are alternative names instead of the main one
Rscript generate/create_fbm.R >generate/create_fbm 2>generate/create_fbm_err # create file based storage
Rscript generate/bootstrap.R >generate/bootstrap_log 2>generate/bootstrap_log_err # calculate the background
Rscript generate/filters.R > generate/filters_log 2> generate/filters_log_err # propagate and save stop words to a file
Rscript generate/package_data.R >generate/package_data_log 2>generate/package_data_log_err # prepare the data to be used, filters, ontologies, caches..

