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

# [x] Output gene-wise contributions to scoring
# [x]--- These numbers may not be very meaningful as they only communicate p-weighted amount of DE.
#        If they could somehow be modified to portray specificity, it would be ideal
# [ ]--- Need to have a condition selector to minimize legend overhead

# [x] Consider interpolating between null distributions

# [ ] Consider result caching
# [ ]--- Make more decisions on what to cache

# [-] Release to lab (ssh -L 12345:localhost:18232 -p 22000 <USERNAME>@willie.msl.ubc.ca)
# [x]--- Think of names

# Known vulnerabilities + weaknesses ----
# 1. Users can use JS to keep sending requests, potentially using all available workers (effectively a very easy DOS attack)
# 2. Users can intentionally crash workers by sending badly formatted URL queries
# 3. Users without JS will have a terrible time
# 4. Users could potentially DOS Gemma by using Gemma DE to query the Gemma API (since requests would be coming from an internal IP address)

DATA.HOLDER$human@gene.meta[gene.Name %in% (strsplit('SORCS3,RBFOX1,GRM5,HIST1H2BN,SHISA9,TCF4,NEGR1,HIST1H3J,DENND1A,DCC,RSRC1,TENM2,TMEM161B,DRD2,PGBD1,ZKSCAN4,HIST1H1B,ERBB4,ZKSCAN8,BTN3A2,PCLO,ZSCAN16,ZSCAN9,TMEM106B,MEF2C,OLFM4,GRM8,ZNF165,LRFN5,OR2B2,HIST1H2BL,ZCCHC7,B3GALTL,BTN2A1,ZSCAN26,RERE,CDH13,ASTN2,CACNA1E,HIST1H2AL,HLA-B,HIST1H4L,ZSCAN12,CHD6,CTNNA3,METTL9,MEGF11,ZSCAN31,ZNF197,KLC1,ZNF660,SPPL3,YLPM1,PCDH9,ZNF445,ZKSCAN7,HIST1H3I,LIN28B,PAX5,PROX2,FAM172A,CELF4,DLST,NRG1,SGCZ,OR12D3,RAB27B,IGSF6,GPC6,PAX6,SGIP1,KDM3A,C16orf45,DCDC5,ESR2,LRP1B,EMILIN3,TRMT61A,XRCC3,SOX5,CNTNAP5,SEMA6D,ANKK1,ZFHX4,LST1,PRSS16,TYR,RFWD2,PQLC2L,BTN1A1,DCDC1,ZDHHC21,TTC12,SDK1,APOPT1,VRK2,CABP1,ZKSCAN3,SAMD5,ADCK3,DENND1B,TAOK3,HS6ST3,MYRF,RTN1,PSORS1C1,CKB,SF3B1,FADS2,GTF2IRD1,NRD1,ZC3H7B,AREL1,RANGAP1,ZNF184,ZDHHC5,HIST1H2BF,FAM120A,KIF15,NKAPL,FCF1,SORBS3,PCDHA2,PCDHA1,PRR34,SCYL1,MR1,BTN3A3,TCTEX1D1,CELF2,CTNND1,HSPA1A,ASCC3,ELAVL2,HIST1H2BO,RPS6KL1,PCDHA3,TRMT10C,ABT1,SCAI,FADS1,CTTNBP2,KMT2A,BEND4,ESRRG,BAZ2B,GPC5,IQCJ-SCHIP1,TCAIM,TMX2,SLC17A3,MED19,ZNF638,CDH22,GRIK5,HARS,HSPE1-MOB4,EP300,HLA-DQB1,PCNP,ZHX3,BCHE,CRB1,C3orf84,MICB,SLC30A9,MARK3,FHIT,MARCH10,CDK14,PLCG1,PSORS1C2,AP3B1,POGZ,TRAF3,CSMD1,TMEM67,PCDHA4,TOPAZ1,PMFBP1,CNTN5,INPP4B,ZNF322,ASIC2,PLA2R1,CHMP3,SOX6,PCDHA5,FANCL,ZNF35,TMEM42,KIAA1143,C11orf31,ACVR1B,ZNF501,RFTN2,TMEM258,TAL1,NICN1,HLA-DQA1,ACTL8,MOB4,CCDC36,PCDHA6,STK19,RHOA,MAP9,FNIP2,RBMS3,PLCL1,SLC44A4,LTBP3,SPRY2,C7orf72,HTT,UBE2M,OTX2,BAG5,CDH9,LPIN3,EPHB2,HMGN4,PPP6C,NOX4,PRR16,EXT1,MGAT4C,EYS,STAU1,HARS2,BAD,MYBPC3,ETFDH,SIM1,FH,ANKS1B,ITPR3,RABEPK,RHOBTB1,BSN,RAB3B,ZNF536,TOP1,CAMKK2,MANEA,ARHGEF25,VPS41,ATP1A3,ITGB6,ASXL3,ANKHD1,PCDHA7,PTPRS,CCS,PHF2,IK,KYNU,PPID,FAM120AOS,ZMAT2,SERPING1,USP3,CACNA2D1,ANKHD1-EIF4EBP3,GINM1,C1QTNF7,MIER1,SLC4A9,PSEN2', ',')[[1]]), entrez.ID] %>% search %>% enrich
