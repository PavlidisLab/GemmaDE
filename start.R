library(shiny)
PROJDIR <- here::here()
DATADIR <- '/cosmos/data/project-data/GemmaDE'
FREEZEDIR <- '/cosmos/data/project-data/GemmaDE/gemma_freeze'
runApp('main', port = 18235, launch.browser = F)

# Roadmap ----
# [?] Consider independent component analysis to reduce feature space
#     per "Content-based microarray search using differential expression profiles"

# [x] Output gene-wise contributions to scoring
# [x]--- These numbers may not be very meaningful as they only communicate p-weighted amount of DE.
#        If they could somehow be modified to portray specificity, it would be ideal
# [ ]--- Need to have a condition selector to minimize legend overhead

# Known vulnerabilities + weaknesses ----
# 1. Users can use JS to keep sending requests, potentially using all available workers (effectively a very easy DOS attack)
# 2. Users can intentionally crash workers by sending badly formatted URL queries
# 3. Users without JS will have a terrible time
# 4. Users could potentially DOS Gemma by using Gemma DE to query the Gemma API (since requests would be coming from an internal IP address)
