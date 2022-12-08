library(MSigDB) #github.com/oganm/MsigDB
source(here::here("main/requirements.R"))
source(here::here("main/dependencies.R"))


tags = getTags(taxa = 'human')

unique_tags = tags %>% dplyr::select(cf.BaseLongUri, cf.ValLongUri) %>% unique
unique_tags %<>% mutate(comparison = glue::glue("{cf.BaseLongUri} vs. {cf.ValLongUri}"))

dnajc5 <- de_search(genes = c('DNAJC5'),
                         taxa = c('human'),geeq = TRUE)


IFNG %>% arrange(desc(`Test Statistic`)) %>% head



MSigDB$human$H %>% lapply(function(x){
  geeq_result =  de_search(genes = x,
                           taxa = c('human'),geeq = TRUE)
  
  no_geeq_result = de_search(genes = x,
                             taxa = c('human'),geeq = FALSE)
  
  list(geeq = geeq_result,
       no_geeq = no_geeq_result)
}) -> results
# saveRDS(results,'msig_db_results.rds')
# results = readRDS('msig_db_results.rds')

# results[[1]]$geeq %>% mutate(Comp = paste0(Baseline, ' vs. ', Value)) %>% {.$Comp == .$`Condition Comparison`} %>% all

MSigDB$human$H %>% names

unique_tags %>% dplyr::filter(grepl("HYPOXIA",toupper(comparison))) %$% comparison %>% unlist %>% unique
unique_tags %>% dplyr::filter(grepl("APOPTOSIS",toupper(comparison))) %$% comparison %>% unlist %>% unique
unique_tags %>% dplyr::filter(grepl("PANCREAS",toupper(comparison))) %$% comparison %>% unlist %>% unique


target_tags = list(HALLMARK_HYPOXIA = c("reference substance role vs. hypoxia",
                          "role vs. hypoxia",
                          "reference subject role vs. hypoxia",
                          "reference subject role vs. 1% O2, hypoxia",
                          "reference subject role vs. 5% O2, hypoxia",
                          "reference subject role vs. Hypoxia",
                          "control vs. hypoxia",
                          "reference substance role vs. 24 h, hypoxia, 0.2% O2",
                          "reference substance role vs. 2 h, hypoxia"),
     
     HALLMARK_APOPTOSIS =c("reference substance role vs. apoptosis inducer",
                           "role vs. apoptosis inducer",
                           "continuant vs. apoptosis inducer",
                           "reference subject role vs. apoptosis inducer"),
     HALLMARK_ESTROGEN_RESPONSE_EARLY = c("reference substance role vs. estrogen",
                                          "role vs. estrogen",
                                          "realizable entity vs. estrogen",
                                          "reference subject role vs. estrogen"),
     HALLMARK_ESTROGEN_RESPONSE_LATE = c("reference substance role vs. estrogen",
                                         "role vs. estrogen",
                                         "realizable entity vs. estrogen",
                                         "reference subject role vs. estrogen"),
     HALLMARK_INTERFERON_ALPHA_RESPONSE = c("reference substance role vs. interferon alpha",
                                            "entity vs. interferon alpha"),
     HALLMARK_INTERFERON_GAMMA_RESPONSE = c("reference substance role vs. interferon gamma"),
     HALLMARK_ALLOGRAFT_REJECTION = unique_tags %>% dplyr::filter(grepl("ALLOGRAFT",toupper(comparison))) %$% comparison %>% unlist %>% unique,
     HALLMARK_ANGIOGENESIS = unique_tags %>% dplyr::filter(grepl("ANGIOGENESIS",toupper(comparison))) %$% comparison %>% unlist %>% unique,
     HALLMARK_COAGULATION = "reference subject role vs. blood coagulation disorder",
     HALLMARK_UV_RESPONSE_UP = c("reference substance role vs. UV-A, radiation exposure",
                                 "reference subject role vs. UV-A",
                                 "reference subject role vs. UV-B",
                                 "reference substance role vs. UV-B, radiation exposure",
                                 "reference subject role vs. Narrow-band UVB, irradiation",
                                 "reference subject role vs. UV radiation, irradiation",
                                 "reference substance role vs. response to UV"),
     HALLMARK_UV_RESPONSE_DN = c("reference substance role vs. UV-A, radiation exposure",
                                 "reference subject role vs. UV-A",
                                 "reference subject role vs. UV-B",
                                 "reference substance role vs. UV-B, radiation exposure",
                                 "reference subject role vs. Narrow-band UVB, irradiation",
                                 "reference subject role vs. UV radiation, irradiation",
                                 "reference substance role vs. response to UV"),
     HALLMARK_INFLAMMATORY_RESPONSE = c("reference subject role vs. inflammation"),
     HALLMARK_ANDROGEN_RESPONSE = c("reference substance role vs. androgen",
                                    "role vs. androgen",
                                    "realizable entity vs. androgen"),
     HALLMARK_PANCREAS_BETA_CELLS = c("exocrine pancreas vs. type B pancreatic cell",
                                      "pancreas vs. type B pancreatic cell"))

MSigDB$human$H %>% names %>% {.[!. %in% names(target_tags)]}


names(target_tags) %>% lapply(function(x){
  geeq = results[[x]]$geeq
  no_geeq = results[[x]]$no_geeq
  
  geeq %<>% arrange(desc(`Test Statistic`))
  no_geeq %<>% arrange(desc(`Test Statistic`))
  
  tags = target_tags[[x]]
  
  geeq_ranks = which(geeq$`Condition Comparison` %in% tags)
  no_geeq_ranks =  which(no_geeq$`Condition Comparison` %in% tags)
  
  geeq_filtered = geeq %>% dplyr::filter(`Condition Comparison` %in% unlist(target_tags))
  no_geeq_filtered = no_geeq %>% dplyr::filter(`Condition Comparison` %in% unlist(target_tags))
  
  geeq_filtered_ranks = which(geeq_filtered$`Condition Comparison` %in% tags)
  no_geeq_filtered_ranks = which(no_geeq_filtered$`Condition Comparison` %in% tags)
  
  list(geeq_ranks = geeq_ranks,
       no_geeq_ranks = no_geeq_ranks,
       geeq_filtered_ranks = geeq_filtered_ranks,
       no_geeq_filtered_ranks = no_geeq_filtered_ranks)
}) -> result_ranks


min_ranks = result_ranks %>% lapply(function(x){sapply(x,min)}) %>% do.call(rbind,.)
mean_ranks = result_ranks %>% lapply(function(x){sapply(x,mean)}) %>% do.call(rbind,.)

wilcox.test(na.omit(mean_ranks[,'geeq_ranks']),na.omit(mean_ranks[,'no_geeq_ranks']),alternative = 'less')
wilcox.test(na.omit(mean_ranks[,'geeq_filtered_ranks']),na.omit(mean_ranks[,'no_geeq_filtered_ranks']),alternative = 'less')

wilcox.test(na.omit(min_ranks[,'geeq_ranks']),na.omit(min_ranks[,'no_geeq_ranks']),alternative = 'less')
wilcox.test(na.omit(min_ranks[,'geeq_filtered_ranks']),na.omit(min_ranks[,'no_geeq_filtered_ranks']),alternative = 'less')

