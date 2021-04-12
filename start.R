source('requirements.R')

source('dependencies.R')
runApp('main', port = 18232, launch.browser = F)

# Roadmap ----
# [?] Consider independent component analysis to reduce feature space
#     per "Content-based microarray search using differential expression profiles"
# [x] Add a biology-related loader

# [x] Support multisessions
# [x]--- Needs to disable the search button
# [x]--- Cancel a process if the client disconnects

# [x] Get a good simulation framework set up for the new analyses
# [x]--- Needs to be tested

# [x] Fix single gene queries
# [x] Fix plot saving
# [x]--- Seems overly complicated and takes awhile
# [x]--- Changed from orca to plotly built-in svg
# [x] Fix table saving

# [/] Output gene-wise contributions to scoring
# [ ]--- These numbers may not be very meaningful as they only communicate p-weighted amount of DE.
#        If they could somehow be modified to portray specificity, it would be ideal
# [ ]--- Need to have a condition selector to minimize legend overhead

# [x] Consider interpolating between null distributions

# [x] Consider result caching
# [ ]--- Make more decisions on what to cache

# [-] Release to lab (ssh -L 12345:localhost:18232 -p 22000 <USERNAME>@willie.msl.ubc.ca)
# [x]--- Think of names

# To be written up ----
# Sex, cell type, tissue specific findings
# Anecdotes (drug-related, metabolic alterations in PD astrocytes)
# Demonstration through simulations
# Breakdown of gaps in knowledge incl. counts for other tools (per species)

# To discuss ----

# 1. BIG poster
# 2. Now have 3x simulations: entirely artificial, spiked in experiment scores and spiked in experiment data
#      Finalizing analysis on spike-in of experiment data
#      Effect of not recalculating prior is non negligible, calculating a few to try to interpolate
# 3. Spent significant time learning about the way these experiments are actually kept in Gemma
# 4. Worked on API wrapper which can soon be unified with the other wrapper (not x-compat but close)
# 5. Dumping more data from API to verify my calculations (some differences)
