# This is a really dirty way to semi-automatically parse the # of data series available on GEO.
# Disclaimer: it doesn't do parsing very intelligently at all, but seems to work for what I need at least

# To use, search on this page: https://www.ncbi.nlm.nih.gov/gds, inputting the following as your search query:
# ("2000/01/01"[Publication Date] : "<YEAR END>/12/31"[Publication Date])
# Replacing <YEAR END> with whatever year you want to query for.
# Under "Study type", show the options "Expression profiling by array" and "Expression profiling by high throughput sequencing"
# Under "Organism", enter as many organisms as you want to quantify (leave them all de-selected but showing)
# Under "Entry type" check "Series" and no other option

# Then check either "Expression profiling by array" or "Expression profiling by high throughput sequencing" (NOT BOTH)
# Press ctrl+A to select the entire page text, ctrl+C to copy it and paste it into GEOpage.txt. Then source this file.

DUMP <- readChar('generate/GEOpage.txt', file.info('generate/GEOpage.txt')$size)

if(!exists('GEO'))
  GEO <- data.table()

GEO <- rbind(GEO,
gsub('.*Search details\n([^\n]*)\nSee more\\.\\.\\.\n.*', '\\1', DUMP) %>% {
  data.table(Year = as.integer(substring(., 24, 27)),
             Type = factor(2 * grepl('high throughput', .) + grepl('array', .), levels = 0:3, labels = c('Nothing', 'Microarray', 'Sequencing', 'Total')) %>%
               as.character)
} %>% cbind(
  strsplit(gsub('\n        ', '~', strsplit(DUMP, 'Organism')[[1]][2]), 'Customize ...', T)[[1]][1] %>%
    strsplit('~') %>% .[[1]] %>% {
      Filter(function(x) x != '', .)
    } %>% lapply(function(x) {
      data.table(Organism = gsub('([^\\(]+)\\(([0-9]+)\\)', '\\1', gsub(',', '', x, fixed = T)),
                 Count = as.integer(gsub('([^\\(]+)\\(([0-9]+)\\)', '\\2', gsub(',', '', x, fixed = T))))
    }) %>% rbindlist %>% rbind(data.table(Organism = 'Total',
                                          Count = as.integer(gsub('.*Series\\(([0-9]+)\\).*', '\\1', gsub(',', '', DUMP, fixed = T)))))
))
