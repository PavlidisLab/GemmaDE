# GEO holdings ----
read.csv(text = 'Year,Quarter,Series,Platforms,Samples
2021,1,"144,666","21,887","4,228,039"
2020,4,"141,740","21,700","4,112,403"
2020,3,"136,638","21,401","3,896,234"
2020,2,"131,878","21,078","3,664,684"
2020,1,"127,183","20,711","3,510,972"
2019,4,"122,753","20,447","3,357,721"
2019,3,"118,590","20,143","3,236,572"
2019,2,"114,770","19,825","3,114,863"
2019,1,"111,015","19,527","2,953,235"
2018,4,"106,808","19,220","2,815,138"
2018,3,"102,879","18,919","2,659,844"
2018,2,"99,695","18,641","2,537,530"
2018,1,"96,434","18,289","2,431,620"
2017,4,"93,033","18,007","2,313,567"
2017,3,"89,405","17,677","2,206,318"
2017,2,"85,890","17,380","2,112,463"
2017,1,"83,054","17,052","2,029,161"
2016,4,"77,671","16,722","1,919,999"
2016,3,"73,555","16,409","1,827,111"
2016,2,"70,696","16,020","1,730,684"
2016,1,"67,104","15,632","1,651,118"
2015,4,"64,076","15,277","1,571,983"
2015,3,"61,175","14,980","1,503,343"
2015,2,"58,626","14,695","1,436,607"
2015,1,"56,020","14,083","1,365,413"
2014,4,"53,622","13,755","1,306,795"
2014,3,"50,985","13,405","1,241,878"
2014,2,"48,534","13,069","1,164,772"
2014,1,"46,334","12,740","1,104,494"
2013,4,"44,100","12,376","1,054,739"
2013,3,"41,756","12,061","1,003,155"
2013,2,"39,234","11,683","949,364"
2013,1,"37,075","11,306","899,088"
2012,4,"34,938","10,889","851,915"
2012,3,"32,848","10,548","807,296"
2012,2,"30,963","10,221","759,246"
2012,1,"29,114","9,925","717,848"
2011,4,"27,342","9,608","673,140"
2011,3,"25,186","9,264","625,389"
2011,2,"23,592","8,836","583,318"
2011,1,"22,050","8,495","541,884"
2010,4,"20,552","8,124","509,535"
2010,3,"18,928","7,811","479,830"
2010,2,"17,494","7,433","448,335"
2010,1,"16,298","7,126","418,410"
2009,4,"14,984","6,752","383,409"
2009,3,"13,740","6,285","353,908"
2009,2,"12,633","5,887","318,072"
2009,1,"11,600","5,633","292,870"
2008,4,"10,671","5,349","272,290"
2008,3,"9,713","4,982","249,468"
2008,2,"8,962","4,734","230,724"
2008,1,"8,224","4,458","210,809"
2007,4,"7,393","4,214","187,389"
2007,3,"6,621","3,795","167,914"
2007,2,"5,944","3,483","151,735"
2007,1,"5,320","3,217","131,412"
2006,4,"4,682","2,876","109,420"
2006,3,"4,228","2,660","99,810"
2006,2,"3,767","2,309","85,270"
2006,1,"3,362","2,134","75,581"
2005,4,"2,869","1,923","63,642"
2005,3,"2,388","1,649","50,228"
2005,2,"2,007","1,394","40,860"
2005,1,"1,763","1,196","33,982"
2004,4,"1,475","1,082","28,258"
2004,3,"1,229",858,"23,116"
2004,2,"1,046",782,"18,347"
2004,1,844,617,"14,936"
2003,4,628,524,"10,386"
2003,3,471,279,"7,498"
2003,2,347,236,"6,049"
2003,1,188,201,"3,972"
2002,4,105,169,"2,645"
2002,3,79,111,"2,103"
2002,2,60,105,"1,874"
2002,1,33,77,"1,187"
2001,4,13,19,670
2001,3,9,14,545
2001,2,2,8,53
2001,1,1,4,48
2000,4,0,2,10
2000,3,0,1,2
') %>% apply(2, function(x) gsub(',', '', x)) %>%
  apply(2, as.numeric) %>% as.data.frame %>%
  mutate(date = paste0(Year, ' Q', Quarter)) -> GEO
GEO %>%
  mutate(YearQuarter = Year + (Quarter - 1) / 4) %>%
  ggplot(aes(YearQuarter, Series)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ poly(x, 4)) +
  theme_classic(base_size = 20) +
  scale_x_continuous(breaks = unique(GEO$Year)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Year', y = '# of Data Series')

# Number of data series exceeds half million in 20 years
GEO %>%
  mutate(YearQuarter = Year + (Quarter - 1) / 4,
         RootSeries = sqrt(Series)) %>% {
           lm(RootSeries ~ YearQuarter, .)
         } %>% predict(data.frame(YearQuarter = 2041)) %>% `^`(2)

CACHE.BACKGROUND$human[distance == 0, .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
CACHE.BACKGROUND$mouse[distance == 0, .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
CACHE.BACKGROUND$rat[distance == 0, .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]

CACHE.BACKGROUND$human[, .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
CACHE.BACKGROUND$mouse[, .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
CACHE.BACKGROUND$rat[, .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]

rbindlist(list(CACHE.BACKGROUND$human,
               CACHE.BACKGROUND$mouse,
               CACHE.BACKGROUND$rat))[distance == 0, .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]


rbindlist(list(CACHE.BACKGROUND$human,
               CACHE.BACKGROUND$mouse,
               CACHE.BACKGROUND$rat))[, .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]

CACHE.BACKGROUND$human

# Overlaps ----
rbindlist(list(CACHE.BACKGROUND$human[distance == 0, .(.N, species = 'Human'),
                                                 .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
               CACHE.BACKGROUND$rat[distance == 0, .(.N, species = 'Rat'),
                                               .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
               CACHE.BACKGROUND$mouse[distance == 0, .(.N, species = 'Mouse'),
                                                 .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)])) %>%
  .[, whichOne := 'Pre-inferrence'] %>%
  .[, mean := mean(N), species] %>%
  .[, median := median(N), species] %>%
  rbind(
    rbindlist(list(CACHE.BACKGROUND$human[, .(.N, species = 'Human'),
                                          .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                   CACHE.BACKGROUND$rat[, .(.N, species = 'Rat'),
                                        .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                   CACHE.BACKGROUND$mouse[, .(.N, species = 'Mouse'),
                                          .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)])) %>%
      .[, whichOne := 'Post-inferrence'] %>%
      .[, mean := mean(N), species] %>%
      .[, median := median(N), species]
  ) %>%
  .[, whichOne := factor(whichOne, levels = c('Pre-inferrence', 'Post-inferrence'))] %>%
  ggplot(aes(log10(N), fill = species)) +
  geom_density() +
  geom_vline(aes(xintercept = log10(mean)), color = 'red') +
  geom_vline(aes(xintercept = log10(median)), color = 'black') +
  facet_grid(vars(whichOne), vars(species), scales = 'free_y', switch = 'y') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = 'Dark2') +
  theme_bw(base_size = 20) +
  theme(legend.position = 'none',
        text = element_text(colour = 'black'), panel.grid = element_blank(),
        strip.placement = 'outside',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = 'bold')) +
  xlab(expression(log[10]~"#"~Overlapping~Experiments)) +
  ylab('Density')

# Fraction of ontology terms ----
rbindlist(list(DATA.HOLDER$human@experiment.meta,
               DATA.HOLDER$mouse@experiment.meta,
               DATA.HOLDER$rat@experiment.meta)) %>% {
                 n <- nrow(.[grepl('http', cf.BaseLongUri) | grepl('http', cf.ValLongUri)])
                 n / nrow(.)
               }

# OE bars ----
search(DATA.HOLDER$human@gene.meta[gene.Name %in% c('XIST', 'RPS4Y1', 'KDM5D', 'EIF1AY', 'DDX3Y'), entrez.ID]) -> searched
enrich(searched) -> enriched

enriched %>% .[, .(id = paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri),
                   Observed = A/B, stat, Empirical = M5, ymin = M5 - S5, ymax = M5 + S5)] %>%
  head(20) %>%
  melt(measure.vars = c('Observed', 'Empirical')) %>%
  setorder(-stat) %>%
  .[, id := factor(id, levels = unique(id), ordered = T)] %>% {
    ggplot(., aes(id, value, fill = variable)) +
      geom_bar(stat = 'identity', position = 'dodge') +
      geom_errorbar(data = copy(.) %>%
                      .[variable == 'Observed', c('ymin', 'ymax') := NA_real_],
                    aes(id, value, group = variable, ymin = ymin, ymax = ymax), width = 0.2, position = position_dodge(0.9)) +
      geom_text(data = copy(.) %>%
                  .[, value := value + 0.0001] %>%
                  .[variable == 'Empirical', value := NA_real_],
                aes(label = round(stat, 1)), fill = NA, position = position_dodge(0.9)) +
      theme_classic(base_size = 20) +
      scale_fill_brewer(palette = 'Dark2', name = 'Group') +
      scale_y_continuous(expand = expansion(add = c(0, 0.001)), labels = scales::percent) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = c(1, 1),
            legend.justification = c(1, 1)) +
      labs(x = element_blank(), y = 'Weighted fraction')
  }
# Sex ----
doSearchEnrich(c('XIST', 'KDM5D', 'RPS4Y1', 'EIF1AY', 'DDX3Y'), 'sex')
rocLike(sex_search_human, sex_enrich_human, truthyCat = quote(cf.Cat == 'biological sex'))
prCurve(sex_search_human, sex_enrich_human, truthyCat = quote(cf.Cat == 'biological sex')) %>% plot

doSearchEnrich(c('A2m', 'Aass', 'Abcb9', 'Abhd3', 'Acot11', 'Acsl6', 'Agt', 'AI464131', 'Akt2', 'Aldh1l1', 'Aldh6a1', 'Aox1', 'Aqp4', 'Aqp9', 'Btbd17', 'Cbs', 'Chrdl1', 'Chst2', 'Cideb', 'Clu', 'Crot', 'Cth', 'Cyp4f14', 'Cyp4f15', 'Cyp7b1', 'Ddhd1', 'Ddo', 'Dhrs11', 'Dio2', 'Egfl6', 'Elovl2', 'Entpd2', 'Eps8', 'Fads2', 'Fgd6', 'Fgfr3', 'Fgfrl1', 'Fjx1', 'Fkbp10', 'Fzd2', 'Fzd7', 'Gas1', 'Gdf10', 'Gfap', 'Gldc', 'Glipr2', 'Gramd3', 'Gstm5', 'Hapln1', 'Hopx', 'Id4', 'Itga7', 'Kctd14', 'Klf15', 'Lfng', 'Lgi4', 'Lix1', 'Mlc1', 'Msmo1', 'Naaa', 'Nat1', 'Nat8', 'Nat8f1', 'Nat8f5', 'Ndp', 'Nr1h3', 'Nrarp', 'Oaf', 'Papss2', 'Paqr7', 'Pax6', 'Pcsk6', 'Pdk4', 'Pdzph1', 'Phkg1', 'Prdx6', 'Pxmp2', 'Retsat', 'Rfx2', 'Rpe65', 'Scara3', 'Sdsl', 'Sema4b', 'Slc12a4', 'Slc14a1', 'Slc22a4', 'Slc27a1', 'Slc7a10', 'Smim3', 'Soat1', 'Sox21', 'St3gal4', 'Tekt1', 'Tmem51', 'Tst', 'Ttpa', 'Ttyh1', 'Ucp3', 'Ugp2', 'Upp1', 'Vit', 'Vnn1', 'Wnt7b', 'Zbtb20', 'Zic5', '1600014C23Rik', '1700003M07Rik', '1700008J07Rik', '1700084J12Rik', '2900052N01Rik', '3110082J24Rik', '4430402I18Rik', '4930515L03Rik', '5930403L14Rik', '6720483E21Rik', '9430041J12Rik', 'A330048O09Rik', 'A330076C08Rik', 'AB124611', 'Acad12', 'Acat3', 'Acsl3', 'Alpk1', 'Apoc3', 'Arhgef26', 'Atp13a4', 'Bcas1os2', 'Bdh2', 'Bglap', 'Bmpr1b', 'C030018K13Rik', 'C030037D09Rik', 'c11_tRNA-Leu-TTG', 'c2_tRNA-Ala-GCG', 'C230035I16Rik', 'Calr4', 'Ccdc190', 'Ccdc69', 'Cfap100', 'Chp2', 'Cml3', 'cx_tRNA-Ala-GCG', 'D630033O11Rik', 'Dmp1', 'Dok7', 'E030003E18Rik', 'Emx2', 'Epo', 'Epx', 'Etnppl', 'Exd1', 'F2rl3', 'Fam181a', 'Fam198a', 'Fam20a', 'Flt3l', 'Gli3', 'Gm10560', 'Gm10658', 'Gm10731', 'Gm12216', 'Gm13031', 'Gm17660', 'Gm5083', 'Gm5475', 'Gm973', 'Gramd1c', 'Hist1h2bq', 'Hist1h2br', 'Hist2h3c1', 'Hk1os', 'Kirrel2', 'Lgr4', 'LOC433374', 'Lonrf3', 'Lpin3', 'Lrrc9', 'Maob', 'Mettl11b', 'Mir99ahg', 'Naprt', 'Nphp3', 'Nr2e1', 'Nwd1', 'Olfr287', 'Olfr288', 'Olfr461', 'Olfr57', 'Olfr688', 'Otx1', 'Paqr8', 'Pax6os1', 'Pex11g', 'Plcd4', 'Plxnb1', 'Ppil6', 'Prr5', 'Psd2', 'Rmst', 'Sctr', 'Sfxn5', 'Shpk', 'Slc25a18', 'Slc39a12', 'Slc6a11', 'Slco1b2', 'Snora69', 'Snord8', 'Sohlh2', 'Stag3', 'Sugct', 'Sycp2', 'Syde2', 'Taf4b', 'Tex11', 'Tmem52', 'Trp63', 'Ugt1a2', 'Ugt1a6a', 'Utp14b', 'Zfp389', 'Zfp783', 'Znrf3', '2810459M11Rik', 'Acaa2', 'Acsbg1', 'Add3', 'Adhfe1', 'Agl', 'Appl2', 'Atp1a2', 'Atp1b2', 'Ccdc80', 'Cd44', 'Cldn10', 'Cyp2d22', 'Cyp2j9', 'Dennd4c', 'Dtna', 'Elovl5', 'Eva1a', 'Ezr', 'F3', 'Fam181b', 'Fmo1', 'Fxyd1', 'Gabrg1', 'Ghr', 'Gjb6', 'Gli2', 'Glud1', 'Gpam', 'Gpld1', 'Grhl1', 'Grin2c', 'Gstm1', 'Hsd11b1', 'Igdcc4', 'Igfbp5', 'Lcat', 'Lhx2', 'Lpp', 'Mfge8', 'Mmp14', 'Nat8f4', 'Ndrg2', 'Ntrk2', 'Ntsr2', 'Oplah', 'Phka1', 'Pla2g7', 'Ppara', 'Ppp1r3c', 'Prelp', 'Prkd1', 'Prodh', 'Ptch1', 'Pygb', 'Rftn2', 'S1pr1', 'Sdc4', 'Serhl', 'Slc15a2', 'Slc1a2', 'Slc1a3', 'Slc1a4', 'Slc38a3', 'Slc41a1', 'Slc4a4', 'Slco1c1', 'Smpd2', 'Smpdl3a', 'Sorcs2', 'Sox9', 'Spag1', 'Tfcp2l1', 'Tlcd1', 'Tmem47', 'Tril', 'Trps1', 'Vav3', 'Vcam1', 'Ccn1'), 'astro', 'mouse')

doSearchEnrich(c('Xist', 'Rps4y1', 'Kdm5d', 'Ddx3y', 'Eif1ay'), 'sex', 'mouse')
rocLike(sex_search_mouse, sex_enrich_mouse, truthyCat = quote(cf.Cat == 'biological sex'), taxa = 'mouse')
prCurve(sex_search_mouse, sex_enrich_mouse, truthyCat = quote(cf.Cat == 'biological sex'), taxa = 'mouse') %>% plot

doSearchEnrich(c('DISP1', 'FGFRL1', 'MSC', 'STRA6', 'TCF21', 'WT1'), 'diaphragm')

doSearchEnrich(c('Sry'), 'sex', 'rat')
rocLike(sex_search_rat, sex_enrich_rat, truthyCat = quote(cf.Cat == 'biological sex'), taxa = 'rat')

# Heart ----
doSearchEnrich(c('MYL7', 'NPPA', 'NPPB'), 'heart')
rocLike(heart_search_human, heart_enrich_human,
        truthyBaseline = quote(grepl('(MYOCD|cardio|cardiac|heart|myocard)', cf.BaseLongUri)),
        truthyValue = quote(grepl('(MYOCD|cardio|cardiac|heart|myocard)', cf.ValLongUri)), max.distance = 0)
prCurve(heart_search_human, heart_enrich_human,
        truthyBaseline = quote(grepl('(MYOCD|cardio|cardiac|heart|myocard|circulator)', cf.BaseLongUri)),
        truthyValue = quote(grepl('(MYOCD|cardio|cardiac|heart|myocard|circulator)', cf.ValLongUri)), max.distance = 0) %>% plot

# Astrocyte ----
doSearchEnrich(c('GFAP', 'SLC1A2', 'ALDH1L1'), 'astro')
doSearchEnrich(c('Gfap'), 'astro', 'mouse')

rocLike(astro_search_human, astro_enrich_human,
        truthyBaseline = quote(grepl('(astro|glia)', cf.BaseLongUri)),
        truthyValue = quote(grepl('(astro|glia)', cf.ValLongUri)), max.distance = 0)
prCurve(astro_search_human, astro_enrich_human,
        truthyBaseline = quote(grepl('(astro)', cf.BaseLongUri)),
        truthyValue = quote(grepl('(astro)', cf.ValLongUri)), max.distance = 0) %>% plot

rocLike(astro_search_mouse, astro_enrich_mouse,
        truthyBaseline = quote(grepl('(astro|glia)', cf.BaseLongUri)),
        truthyValue = quote(grepl('(astro|glia)', cf.ValLongUri)), taxa = 'mouse', max.distance = 0)
prCurve(astro_search_mouse, astro_enrich_mouse,
        truthyBaseline = quote(grepl('(astro|glia)', cf.BaseLongUri)),
        truthyValue = quote(grepl('(astro|glia)', cf.ValLongUri)), taxa = 'mouse', max.distance = 0) %>% plot

# Oligo ----
doSearchEnrich(c('PDGFA', 'CSPG4', 'OLIG1', 'OLIG2', 'OLIG3'), 'oligo')
rocLike(oligo_search_human, oligo_enrich_human,
        truthyBaseline = quote(grepl('(oligo|glia)', cf.BaseLongUri)),
        truthyValue = quote(grepl('(oligo|glia)', cf.ValLongUri)), max.distance = 0)

# Drugs ----
doSearchEnrich(c('LHCGR', 'GNRHR'), 'goserelin')
rocLike(goserelin_search_human, goserelin_enrich_human,
        truthyBaseline = quote(grepl('goserelin', cf.BaseLongUri, ignore.case = T)),
        truthyValue = quote(grepl('goserelin', cf.ValLongUri, ignore.case = T)))

# Helpers ----
doSearchEnrich <- function(genes, name, taxa = getConfig('taxa')$value) {
  tmp <- search(DATA.HOLDER[[taxa]]@gene.meta[gene.Name %in% genes, entrez.ID], getConfig(taxa = taxa))
  assign(paste0(name, '_search_', taxa), tmp, envir = globalenv())
  
  tmp <- enrich(tmp, getConfig(taxa = taxa))
  assign(paste0(name, '_enrich_', taxa), tmp, envir = globalenv())
}

rocLike <- function(search, enrich, truthyBaseline = quote(F),
                    truthyValue = truthyBaseline, truthyCat = quote(F),
                    taxa = getConfig('taxa')$value, max.distance = Inf) {
  search %>%
    merge(DATA.HOLDER[[taxa]]@experiment.meta[, .(cf.Cat, cf.BaseLongUri = cf.Baseline,
                                                  cf.ValLongUri = cf.Val, rsc.ID, ee.ID, n.DE, ee.qScore)],
          by.x = 'rn', by.y = 'rsc.ID', sort = F) %>%
    .[, I := .I / max(.I)] %>%
    .[, isGood := eval(truthyBaseline) | eval(truthyValue) | eval(truthyCat)] %>%
    .[, isGood := cumsum(isGood) / sum(isGood)] %>%
    .[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, I, isGood, type = 'Search')] %>% {
      print(head(., 50))
      .
    } %>%
    rbind(enrich[distance <= max.distance, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
            .[, I := .I / max(.I)] %>%
            .[, isGood := eval(truthyBaseline) | eval(truthyValue) | eval(truthyCat)] %>%
            .[, isGood := cumsum(isGood) / sum(isGood)] %>%
            .[, type := 'Enrich']) %>%
    ggplot(aes(I, isGood, color = type)) + geom_line(size = 2) +
    theme_classic(base_size = 20) + scale_color_brewer(palette = 'Dark2', name = 'Module') +
    labs(x = 'Fractional index', y = 'Cumulative distribution') +
    geom_abline(slope = 1, intercept = 0, size = 0.5) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0)), labels = scales::percent) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), labels = scales::percent)
}

prCurve <- function(search, enrich, truthyBaseline = quote(F),
                    truthyValue = truthyBaseline, truthyCat = quote(F),
                    taxa = getConfig('taxa')$value, max.distance = Inf, whichOne = 'Enrich') {
  search %>%
    merge(DATA.HOLDER[[taxa]]@experiment.meta[, .(cf.Cat, cf.BaseLongUri = cf.Baseline,
                                                  cf.ValLongUri = cf.Val, rsc.ID, ee.ID, n.DE, ee.qScore)],
          by.x = 'rn', by.y = 'rsc.ID', sort = F) %>%
    .[, I := .I / max(.I)] %>%
    .[, isGood := eval(truthyBaseline) | eval(truthyValue) | eval(truthyCat)] %>%
    .[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, I, isGood, score, type = 'Search')] %>% {
      print(head(., 50))
      .
    } %>%
    rbind(enrich[distance <= max.distance, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, stat)] %>%
            .[, I := .I / max(.I)] %>%
            .[, isGood := eval(truthyBaseline) | eval(truthyValue) | eval(truthyCat)] %>%
            .[, score := stat] %>%
            .[, stat := NULL] %>%
            .[, type := 'Enrich']) %>% {
              print(.[, sum(isGood), type])
              pr.curve(scores.class0 = .[type == whichOne & isGood == T, score],
                       scores.class1 = .[type == whichOne & isGood == F, score], curve = T)
            }
}
