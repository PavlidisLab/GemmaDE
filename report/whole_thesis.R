# this is an attempt to have an up to date working version of all plots and claims
# included in Jordan's thesis, in order of appearance -Ogan

devtools::load_all()
source(here::here("main/requirements.R"))
source(here::here("main/dependencies.R"))

# Chapter 1, figure 1.2 GEO and Gemma Contents----------
experiments = readRDS('report/experiments.rds')

# did not regenerate this file from scratch. generate/quantifyGEO.R deals with it
# and needs to be tested
GEO <- fread('generate//GEO_holdings_2021.csv')
GEO %>% melt(measure.vars = 'Count') %>%
  rbind(GEO %>% melt(measure.vars = 'Count') %>%
          .[Organism %in% c('Homo sapiens', 'Mus musculus', 'Rattus norvegicus'), .(value = sum(value)), .(Year, Type)] %>% {
            data.table(Year = .$Year,
                       Type = .$Type,
                       Organism = 'Mammalian',
                       variable = 'Count',
                       value = .$value)
          }) %>%
  .[, .(value = sum(value)), .(Year, Organism)] %>% {
    data.table(Year = .$Year,
               Type = 'Aggregate',
               Organism = .$Organism,
               variable = 'Count',
               value = .$value)
  } %>%
  rbind(GEO %>% melt(measure.vars = 'Count') %>%
          rbind(GEO %>% melt(measure.vars = 'Count') %>%
                  .[Organism %in% c('Homo sapiens', 'Mus musculus', 'Rattus norvegicus'), .(value = sum(value)), .(Year, Type)] %>% {
                    data.table(Year = .$Year,
                               Type = .$Type,
                               Organism = 'Mammalian',
                               variable = 'Count',
                               value = .$value)
                  })) %>%
  #.[Organism == 'Mammalian'] %>%
  #.[, ingemma := experiments[taxon.Name %in% c('human', 'mouse', 'rat'), length(unique(ee.ID))]] %>%
  .[, Type := factor(Type, levels = c('Microarray', 'Sequencing', 'Aggregate'), labels = c('Microarray', 'Sequencing', 'Total'), ordered = T)] %>%
  .[Organism %in% c('Homo sapiens', 'Mus musculus', 'Rattus norvegicus')] %>%
  .[, ingemma := case_when(Organism == 'Homo sapiens' ~ experiments[experiment.Database == 'GEO' & taxon.Name == 'human' &
                                                                      (Type == 'Total' | (platform.Type %in% c('GENELIST', 'SEQUENCING')) == (Type == 'Sequencing')), length(unique(experiment.ID))],
                           Organism == 'Mus musculus' ~ experiments[experiment.Database == 'GEO' & taxon.Name == 'mouse' &
                                                                      (Type == 'Total' | (platform.Type %in% c('GENELIST', 'SEQUENCING')) == (Type == 'Sequencing')), length(unique(experiment.ID))],
                           Organism == 'Rattus norvegicus' ~ experiments[experiment.Database == 'GEO' & taxon.Name == 'rat' &
                                                                           (Type == 'Total' | (platform.Type %in% c('GENELIST', 'SEQUENCING')) == (Type == 'Sequencing')), length(unique(experiment.ID))]), Type] %>%
  ggplot(aes(Year, value, color = Type)) +
  geom_point() + geom_line(lwd = 1.1) +
  #geom_hline(aes(yintercept = ingemma), color = 'black', linetype = 'dashed', lwd = 1.5) +
  geom_hline(aes(color = Type, yintercept = ingemma), lwd = 1.4, linetype = 'dashed') +
  facet_wrap(~Organism, scales = 'free_y', ncol = 1) +
  theme_bw(20) + xlab(element_blank()) + ylab('# of Data Series') +
  theme(strip.background = element_blank(), strip.text = element_text(face = 'bold.italic', hjust = 0)) +
  scale_color_brewer(palette = 'Dark2') +
  scale_x_continuous(expand = c(0, 0))


# Chapter 2 Figure 2.2 --------
rbindlist(lapply(names(CACHE.BACKGROUND), function(i) data.table(taxon = i, CACHE.BACKGROUND[[i]]))) %>%
  .[, rsc.ID := as.integer(as.factor(rsc.ID))] -> mCache

subjects = list(
  System = c('nervous', 'reproductive', 'digestive', 'respiratory', 'hemolymphoid',
             'endocrine', 'exocrine', 'cardiovascular', 'hepatobilliary'),
  `Organ/Tissue` = c('spinal cord','kidney','spleen','heart','intestine','muscle',
                     'lung','liver','blood','brain'),
  `Cell Type` = c(
    'microglial cell','embryonic stem cell','B cell','T cell','macrophage',
    'fibroblast','glial cell','epithelial cell','leukocyte','hematopoietic cell'
  )
)


lapply(subjects %>% unlist,
       function(topic) {
         data.table(Topic = topic, mCache[grepl(topic, cf.BaseLongUri, T) | grepl(topic, cf.ValLongUri, T), .(EE = length(unique(ee.ID)), RSC = length(unique(rsc.ID))), taxon])
       }) %>% rbindlist -> tagSummary


tagSummary %>%
  .[, subject := tagSummary$Topic %>% sapply(function(x){subjects %>% sapply(function(y){x %in% y}) %>% which %>% names})
  ] %>%
  melt(measure.vars = c('EE', 'RSC')) %>%
  .[!is.na(value) & variable == 'RSC'] %>%
  .[, stat := sum(value), Topic] %>%
  setorder(-stat) %>%
  .[, Topic := factor(Topic, levels = unique(Topic), ordered = T)] %>%
  ggplot(aes(Topic, value, fill = taxon)) +
  geom_bar(stat = 'identity') + facet_wrap(~subject, scales = 'free') + coord_flip() +
  theme_classic(20) + scale_y_continuous(expand = c(0, 0)) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = 'bold'), legend.text.align = 0, legend.justification = 'top') +
  xlab(element_blank()) + ylab('Condition Comparisons') +
  scale_fill_brewer(palette = 'Dark2', name = 'Taxon', labels = c(expression(italic('H. sapiens')), expression(italic('M. musculus')), expression(italic('R. norvegicus'))))


# Chapter 3 Figure 3.1 Overlaps -----------
rbindlist(list(DATA.HOLDER$human@experiment.meta[, .(N = length(unique(ee.ID)), species = 'Human'),
                                                 .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
               DATA.HOLDER$rat@experiment.meta[, .(N = length(unique(ee.ID)), species = 'Rat'),
                                               .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
               DATA.HOLDER$mouse@experiment.meta[, .(N = length(unique(ee.ID)), species = 'Mouse'),
                                                 .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)])) %>%
  .[, whichOne := 'Pre-inference'] %>%
  .[, mean := mean(N), species] %>%
  .[, median := median(N), species] %>%
  rbind(
    rbindlist(list(CACHE.BACKGROUND$human[, .(N = length(unique(ee.ID)), species = 'Human'),
                                          .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                   CACHE.BACKGROUND$rat[, .(N = length(unique(ee.ID)), species = 'Rat'),
                                        .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                   CACHE.BACKGROUND$mouse[, .(N = length(unique(ee.ID)), species = 'Mouse'),
                                          .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)])) %>%
      .[, whichOne := 'Post-inference'] %>%
      .[, mean := mean(N), species] %>%
      .[, median := median(N), species]
  ) %>%
  .[, whichOne := factor(whichOne, levels = c('Pre-inference', 'Post-inference'))] %>% {
    print(.[, .(whichOne, mean, median)] %>% unique)
    .
  } %>%
  ggplot(aes(log10(N), fill = species)) +
  geom_density() +
  geom_vline(aes(xintercept = log10(mean)), color = 'red') +
  geom_vline(aes(xintercept = log10(median)), color = 'black') +
  facet_grid(vars(whichOne), vars(species), scales = 'free_y', switch = 'y') +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2.25), breaks = c(0, log10(c(2, 5)), 1, log10(c(25, 50)), 2), labels = function(x) round(10^x)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = 'Dark2') +
  theme_bw(base_size = 20) +
  theme(legend.position = 'none',
        text = element_text(colour = 'black'), panel.grid = element_blank(),
        strip.placement = 'outside',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = 'bold')) +
  xlab(expression("#"~Overlapping~Experiments)) +
  ylab('Density')


# chapter 3.2 demonstration of gene list results -----------
sex_results <- de_search(genes = c('RPS4Y1', 'EIF1AY', 'DDX3Y', 'KDM5D', 'XIST'),
                         taxa = c('human'),geeq = TRUE)

sex_results_nogeeq <- de_search(genes = c('RPS4Y1', 'EIF1AY', 'DDX3Y', 'KDM5D', 'XIST'),
                                taxa = c('human'),geeq = FALSE)

ct_markers <- c("Cyp4f15", "Grin2c", "Col4a5", "Heph", "Celsr1", "Egfr", "Lgi4", 
                "Slc38a3", "Aass", "Fkbp10", "Slc7a10", "Cyp2d22", "Cd38", "Cyp4f14", 
                "Cbs", "Slc1a2", "Emp2", "Axl", "Slc14a1", "Btbd17")

cell_type_results <- de_search(genes = ct_markers,
                               taxa = 'mouse',geeq = TRUE)

cell_type_results_nogeeq <- de_search(genes = ct_markers,
                                      taxa = 'mouse',geeq = FALSE)

glucocorticoid_genes <- c("BCL6", "BIRC3", "CEBPD", "ERRFI1", "FBXL16", "FKBP5", "GADD45B", 
                          "IRS2", "KLF9", "PDK4", "PER1", "RGCC", "RGS2", "SEC14L2", "SLC16A12", 
                          "TFCP2L1", "TSC22D3")

# metadata seems to have changed a bit splitting left ventricular assist device from ventricular assist device in analysis. these are plain text annotations
CACHE.BACKGROUND$human$cf.ValLongUri[grepl('assist device',CACHE.BACKGROUND$human$cf.ValLongUri)] = 'ventricular assist device'

glucocorticoid_results <- de_search(genes = glucocorticoid_genes,
                                    taxa = 'human',geeq = TRUE)

glucocorticoid_results_nogeeq <- de_search(genes = glucocorticoid_genes,
                                           taxa = 'human',geeq = FALSE)


heart_results <- de_search(genes = c('MYL7', 'NPPA', 'NPPB'),taxa = 'human',geeq=TRUE)
heart_results_nogeeq <- de_search(genes = c('MYL7', 'NPPA', 'NPPB'),taxa = 'human',geeq= FALSE)



scz_genes = c("HARS", "HARS2", "HCN1", "HIRIP3", "HSPA9", "HSPD1", "HSPE1", 
              "IGSF9B", "IK", "IMMP2L", "INA", "INO80E", "IREB2", "ITIH1", 
              "ITIH3", "ITIH4", "KCNB1", "KCNJ13", "KCNV1", "KCTD13", "KDM3B", 
              "KDM4A", "KLC1", "L3MBTL2", "LCAT", "LRP1", "LRRIQ3", "LUZP2", 
              "MAD1L1", "MAN2A1", "MAN2A2", "MAPK3", "MARS2", "MAU2", "MDK", 
              "MED19", "MEF2C", "MIR137", "MIR548AJ2", "MMP16", "MPHOSPH9", 
              "MPP6", "MSANTD2", "MSL2", "MUSTN1", "MYO15A", "MYO1A", "NAB2", 
              "NAGA", "NCAN", "NCK1", "NDUFA13", "NDUFA2", "NDUFA4L2", "NDUFA6", 
              "NEK1", "NEK4", "NFATC3", "NGEF", "NISCH", "NLGN4X", "NOSIP", 
              "NRGN", "NRN1L", "NT5C2", "NT5DC2", "NUTF2", "NXPH4", "OGFOD2", 
              "OSBPL3", "OTUD7B", "PAK6", "PARD6A", "PBRM1", "PBX4", "PCCB", 
              "PCDHA1", "PCDHA10", "PCDHA2", "PCDHA3", "PCDHA4", "PCDHA5", 
              "PCDHA6", "PCDHA7", "PCDHA8", "PCDHA9", "PCGEM1", "PCGF6", "PDCD11", 
              "PITPNM2", "PJA1", "PLA2G15", "PLCB2", "PLCH2", "PLCL1", "PLEKHO1", 
              "PODXL", "PPP1R13B", "PPP1R16B", "PPP2R3A", "PPP4C", "PRKD1", 
              "PRR12", "PRRG2", "PSKH1", "PSMA4", "PSMB10", "PSMD6", "PTGIS", 
              "PTN", "PTPRF", "PUS7", "R3HDM2", "RAI1", "RANBP10", "RANGAP1", 
              "RCN3", "REEP2", "RERE", "RFTN2", "RGS6", "RILPL2", "RIMS1", 
              "RRAS", "SATB2", "SBNO1", "SCAF1", "SDCCAG8", "SERPING1", "SEZ6L2", 
              "SF3B1", "SFXN2", "SGSM2", "SHISA8", "SHMT2", "SLC12A4", "SLC32A1", 
              "SLC35G2", "SLC38A7", "SLC39A8", "SLC45A1", "SLC4A10", "SLC7A6", 
              "SLC7A6OS", "SMDT1", "SMG6", "SMIM4", "SNAP91", "SNX19", "SPCS1", 
              "SREBF1", "SREBF2", "SRPK2", "SRR", "STAB1", "STAC3", "STAG1", 
              "STAT6", "SUGP1", "TAC3", "TAF5", "TAOK2", "TBC1D5", "TBX6", 
              "TCF20", "TCF4", "THAP11", "THOC7", "TLE1", "TLE3", "TM6SF2", 
              "TMCO6", "TMEM110-MUSTN1", "TMEM219", "TMTC1", "TMX2", "TNFRSF13C", 
              "TOM1L2", "TRANK1", "TRIM8", "TRMT61A,ABCB9", "ACD", "ACTR5", 
              "ADAMTSL3", "AKT3", "ALDOA", "AMBRA1", "ANKRD44", "ANKRD63", 
              "ANP32E", "APH1A", "APOPT1", "ARHGAP1", "ARL3", "ARL6IP4", "AS3MT", 
              "ASPHD1", "ATG13", "ATP2A2", "ATPAF2", "ATXN7", "BAG5", "BCL11B", 
              "BOLL", "BTBD18", "C11orf87", "C12orf42", "C12orf65", "C16orf86", 
              "C16orf92", "C1orf54", "C2orf69", "C3orf49", "CA14", "CA8", "CACNA1C", 
              "CACNA1I", "CACNB2", "CCDC39", "CD14", "CD46", "CDC25C", "CDK2AP1", 
              "CENPM", "CENPT", "CHADL", "CHRM4", "CHRNA3", "CHRNA5", "CHRNB4", 
              "CILP2", "CKAP5", "CKB", "CLCN3", "CLP1", "CLU", "CNKSR2", "CNNM2", 
              "CNOT1", "CNTN4", "COQ10B", "CR1L", "CREB3L1", "CSMD1", "CTNNA1", 
              "CTNND1", "CTRL", "CUL3", "CYP17A1", "CYP26B1", "CYP2D6", "DDX28", 
              "DGKI", "DGKZ", "DNAJC19", "DND1", "DOC2A", "DPEP2", "DPEP3", 
              "DPP4", "DPYD", "DRD2", "DRG2", "EDC4", "EFHD1", "EGR1", "ENKD1", 
              "EP300", "EPC2", "EPHX2", "ERCC4", "ESAM", "ESRP2", "ETF1", "F2", 
              "FAM53C", "FAM57B", "FANCL", "FES", "FURIN", "FUT9", "FXR1", 
              "GALNT10", "GATAD2A", "GDPD3", "GFOD2", "GFRA3", "GID4", "GIGYF2", 
              "GLT8D1", "GNL3", "GOLGA6L4", "GPM6A", "GRAMD1B", "GRIA1", "GRIN2A", 
              "GRM3", "HAPLN4", "HARBI1")

scz_results <- de_search(genes = scz_genes,
                         taxa = 'human', geeq = TRUE)

scz_results_nogeeq <- de_search(genes = scz_genes,
                                taxa = 'human', geeq = FALSE)

depression_genes = c("SORCS3", "RBFOX1", "GRM5", "HIST1H2BN", "SHISA9", "TCF4", 
                     "NEGR1", "HIST1H3J", "DENND1A", "DCC", "RSRC1", "TENM2", "TMEM161B", 
                     "DRD2", "PGBD1", "ZKSCAN4", "HIST1H1B", "ERBB4", "ZKSCAN8", "BTN3A2", 
                     "PCLO", "ZSCAN16", "ZSCAN9", "TMEM106B", "MEF2C", "OLFM4", "GRM8", 
                     "ZNF165", "LRFN5", "OR2B2", "HIST1H2BL", "ZCCHC7", "B3GALTL", 
                     "BTN2A1", "ZSCAN26", "RERE", "CDH13", "ASTN2", "CACNA1E", "HIST1H2AL", 
                     "HLA-B", "HIST1H4L", "ZSCAN12", "CHD6", "CTNNA3", "METTL9", "MEGF11", 
                     "ZSCAN31", "ZNF197", "KLC1", "ZNF660", "SPPL3", "YLPM1", "PCDH9", 
                     "ZNF445", "ZKSCAN7", "HIST1H3I", "LIN28B", "PAX5", "PROX2", "FAM172A", 
                     "CELF4", "DLST", "NRG1", "SGCZ", "OR12D3", "RAB27B", "IGSF6", 
                     "GPC6", "PAX6", "SGIP1", "KDM3A", "C16orf45", "DCDC5", "ESR2", 
                     "LRP1B", "EMILIN3", "TRMT61A", "XRCC3", "SOX5", "CNTNAP5", "SEMA6D", 
                     "ANKK1", "ZFHX4", "LST1", "PRSS16", "TYR", "RFWD2", "PQLC2L", 
                     "BTN1A1", "DCDC1", "ZDHHC21", "TTC12", "SDK1", "APOPT1", "VRK2", 
                     "CABP1", "ZKSCAN3", "SAMD5", "ADCK3", "DENND1B", "TAOK3", "HS6ST3", 
                     "MYRF", "RTN1", "PSORS1C1", "CKB", "SF3B1", "FADS2", "GTF2IRD1", 
                     "NRD1", "ZC3H7B", "AREL1", "RANGAP1", "ZNF184", "ZDHHC5", "HIST1H2BF", 
                     "FAM120A", "KIF15", "NKAPL", "FCF1", "SORBS3", "PCDHA2", "PCDHA1", 
                     "PRR34", "SCYL1", "MR1", "BTN3A3", "TCTEX1D1", "CELF2", "CTNND1", 
                     "HSPA1A", "ASCC3", "ELAVL2", "HIST1H2BO", "RPS6KL1", "PCDHA3", 
                     "TRMT10C", "ABT1", "SCAI", "FADS1", "CTTNBP2", "KMT2A", "BEND4", 
                     "ESRRG", "BAZ2B", "GPC5", "IQCJ-SCHIP1", "TCAIM", "TMX2", "SLC17A3", 
                     "MED19", "ZNF638", "CDH22", "GRIK5", "HARS", "HSPE1-MOB4", "EP300", 
                     "HLA-DQB1", "PCNP", "ZHX3", "BCHE", "CRB1", "C3orf84", "MICB", 
                     "SLC30A9", "MARK3", "FHIT", "MARCH10", "CDK14", "PLCG1", "PSORS1C2", 
                     "AP3B1", "POGZ", "TRAF3", "CSMD1", "TMEM67", "PCDHA4", "TOPAZ1", 
                     "PMFBP1", "CNTN5", "INPP4B", "ZNF322", "ASIC2", "PLA2R1", "CHMP3", 
                     "SOX6", "PCDHA5", "FANCL", "ZNF35", "TMEM42", "KIAA1143", "C11orf31", 
                     "ACVR1B", "ZNF501", "RFTN2", "TMEM258", "TAL1", "NICN1", "HLA-DQA1", 
                     "ACTL8", "MOB4", "CCDC36", "PCDHA6", "STK19", "RHOA", "MAP9", 
                     "FNIP2", "RBMS3", "PLCL1", "SLC44A4", "LTBP3", "SPRY2", "C7orf72", 
                     "HTT", "UBE2M", "OTX2", "BAG5", "CDH9", "LPIN3", "EPHB2", "HMGN4", 
                     "PPP6C", "NOX4", "PRR16", "EXT1", "MGAT4C", "EYS", "STAU1", "HARS2", 
                     "BAD", "MYBPC3", "ETFDH", "SIM1", "FH", "ANKS1B", "ITPR3", "RABEPK", 
                     "RHOBTB1", "BSN", "RAB3B", "ZNF536", "TOP1", "CAMKK2", "MANEA", 
                     "ARHGEF25", "VPS41", "ATP1A3", "ITGB6", "ASXL3", "ANKHD1", "PCDHA7", 
                     "PTPRS", "CCS", "PHF2", "IK", "KYNU", "PPID", "FAM120AOS", "ZMAT2", 
                     "SERPING1", "USP3", "CACNA2D1", "ANKHD1-EIF4EBP3", "GINM1", "C1QTNF7", 
                     "MIER1", "SLC4A9", "PSEN2")

depression_results <- de_search(genes = depression_genes,
                                taxa = 'human',geeq = TRUE)
depression_results %<>% arrange(desc(`Test Statistic`))
depression_results_nogeeq <- de_search(genes = depression_genes,
                                       taxa = 'human',geeq = FALSE)
depression_results_nogeeq %<>% arrange(desc(`Test Statistic`))

# 3.2.1 Sex specific genes -----

# first hit is male vs. female and least contribution is RPS4Y1
sex_results %>% dplyr::arrange(desc(`Test Statistic`)) %>% {.[1,]}

# RPS4Y1 alone has the strongest link on its own
RPS4Y1 <- de_search(genes = c('RPS4Y1'),
                    taxa = c('human'),geeq = TRUE)

EIF1AY <- de_search(genes = c('EIF1AY'),
                    taxa = c('human'),geeq = TRUE)

# 3.2.2 Cell type specific genes ------
# top 20 astrocyte markers ranked by their main silhouette coefficient
# not needed but returns different genes than jordan's. not sure what the discrepancy is
# results are similar
# astro_genes <- readr::read_delim(delim = ' ',
#                                  col_names = c('Gene','fc','mainsil','minsil'),
#                                  'https://raw.githubusercontent.com/PavlidisLab/neuroExpressoAnalysis/master/analysis/01.SelectGenes/Markers_1.2/Quick/All_CellTypes/Astrocyte')  %>%
#   arrange(desc(mainsil)) %$% Gene %>% {.[1:20]}
# 
# astro_results <- de_search(ct_markers,taxa = 'mouse',geeq = TRUE)
# astro_results_no_geeq <- de_search(ct_markers,taxa = 'mouse',geeq = FALSE)


cell_type_results %>% dplyr::arrange(desc(`Test Statistic`)) %>% {.[1:50,]} %>% 
  dplyr::mutate(`Condition Comparison` = factor(`Condition Comparison`,levels = `Condition Comparison`)) %>% 
  ggplot(aes(y = `Condition Comparison`, x = `Test Statistic`,fill = `Ontology Steps`)) + geom_bar(stat = 'identity')

# 3.2.3 drug modulated genes ----------
# this result no longer appears to replicate
glucocorticoid_results


# 3.2.5.1 Gemma DE can offer complementary insights to GO enrichment analysis -----

# a small addition to see the effect of geeq on the results here
scz_results_no_geeq <- de_search(genes = scz_genes,
                                 taxa = 'human', geeq = FALSE)

library(ermineR)

go_enrichment = ermineR::ora(hitlist = scz_genes,annotation = 'Generic_human')
View(go_enrichment$results)

