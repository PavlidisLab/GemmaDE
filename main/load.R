fixOntoGenes <- function() {
  lapply(getConfig(key = "taxa")$core, function(taxa) {
    matches <- do.call(
      rbind,
      str_match_all(
        DATA.HOLDER[[taxa]]@experiment.meta[, c(
          as.character(cf.BaseLongUri),
          as.character(cf.ValLongUri)
        )],
        "http://purl.org/commons/record/ncbi_gene/(\\d*)"
      )
    ) %>% unique()

    matches[, 2] <- lapply(getConfig(key = "taxa")$core, function(taxa) {
      DATA.HOLDER[[taxa]]@gene.meta[
        entrez.ID %in% matches[, 2],
        .(entrez.ID, gene.Name = paste0(gene.Name, " [", taxa, "]"))
      ]
    }) %>%
      data.table::rbindlist() %>%
      unique(by = "entrez.ID") %>%
      .[match(matches[, 2], entrez.ID), gene.Name]

    data.table::data.table(Node_Long = matches[, 1], Definition = matches[, 2], OntologyScope = "TGEMO")
  }) %>%
    data.table::rbindlist() %>%
    {
      rbind(ONTOLOGIES.DEFS[, .(
        Node_Long = as.character(Node_Long),
        Definition = as.character(Definition),
        OntologyScope = as.character(OntologyScope)
      )], .) %>%
        .[, c("Node_Long", "Definition", "OntologyScope") := list(
          as.factor(Node_Long),
          as.factor(Definition),
          as.factor(OntologyScope)
        )]
    }
}


setClass("EData", representation(
  taxon = "character", data = "list",
  experiment.meta = "data.table", gene.meta = "data.table",
  go = "data.table"
))

# Load the data into the global environment
if (!exists("ONTOLOGIES") || !exists("ONTOLOGIES.DEFS")) {
  ONTOLOGIES <- data.table::fread("/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED.TSV")
  ONTOLOGIES.DEFS <- data.table::fread("/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED_DEF.TSV")

  ONTOLOGIES[, c("ChildNode", "ParentNode")] <- NULL
  ONTOLOGIES.DEFS$Node <- NULL

  ONTOLOGIES$ChildNode_Long <- ONTOLOGIES$ChildNode_Long %>% as.factor()
  ONTOLOGIES$ParentNode_Long <- ONTOLOGIES$ParentNode_Long %>% as.factor()
  ONTOLOGIES$RelationType <- ONTOLOGIES$RelationType %>% as.factor()
  ONTOLOGIES$OntologyScope <- ONTOLOGIES$OntologyScope %>% as.factor()

  ONTOLOGIES.DEFS$Node_Long <- ONTOLOGIES.DEFS$Node_Long %>% as.factor()
  ONTOLOGIES.DEFS$OntologyScope <- ONTOLOGIES.DEFS$OntologyScope %>% as.factor()
}

.DATA_PATH <- paste(DATADIR, 'DATA.HOLDER.light.rds', sep='/')

# Load the lite versions if they're already created.
if (!exists("DATA.HOLDER")) {
  if (file.exists(.DATA_PATH)) {
    DATA.HOLDER <- readRDS(.DATA_PATH)
  } else {
    stop("Couldn't find DATA.HOLDER, run generate/fbm.R first.")
  }

  DATA.HOLDER$artificial <- NULL

  ONTOLOGIES.DEFS <- fixOntoGenes()
}

# Read existing FBMs
for (i in names(DATA.HOLDER)) {
  DATA.HOLDER[[i]]@data$zscore <- bigstatsr::big_attach(paste0(paste(DATADIR, 'fbm/', sep='/'), i, "/zscores.rds"))
  attr(DATA.HOLDER[[i]]@data$zscore, ".dimnames") <- readRDS(paste0(paste(DATADIR, 'fbm/', sep='/'), i, "/z.dimnames.rds"))

  DATA.HOLDER[[i]]@data$adj.pv <- bigstatsr::big_attach(paste0(paste(DATADIR, 'fbm/', sep='/'), i, "/adjpvs.rds"))
  attr(DATA.HOLDER[[i]]@data$adj.pv, ".dimnames") <- readRDS(paste0(paste(DATADIR, 'fbm/', sep='/'), i, "/p.dimnames.rds"))
}
rm(i)

gc()

dimnames.FBM <- function(object, ...) {
  attr(object, ".dimnames")
}

if (!exists("CACHE.BACKGROUND")) {
  # Pre-load all ontology expansions
  if (file.exists(paste(DATADIR, 'CACHE.BACKGROUND.rds', sep='/'))) {
    CACHE.BACKGROUND <- readRDS(paste(DATADIR, 'CACHE.BACKGROUND.rds', sep='/'))
  } else {
    CACHE.BACKGROUND <- lapply(names(DATA.HOLDER), precomputeTags) %>%
      setNames(names(DATA.HOLDER))

    saveRDS(CACHE.BACKGROUND, paste(DATADIR, 'CACHE.BACKGROUND.rds', sep='/'))
  }
}

if (!exists("NULLS")) {
  NULLS <- lapply(names(DATA.HOLDER), function(taxa) {
    tryCatch(readRDS(paste0(paste(DATADIR, 'updated_nulls/', sep='/'), taxa, ".rds")) %>%
      .[, .(rn, score.mean, score.sd)], error = function(e) NULL)
  }) %>%
    `names<-`(names(DATA.HOLDER)) %>%
    {
      Filter(Negate(is.null), .)
    }
}

# Compile all unique gene names to use as choices in the search bar
ALL.GENES <- list()
for (taxon in names(DATA.HOLDER)) {
  ALL.GENES[[taxon]] <- DATA.HOLDER[[taxon]]@gene.meta$gene.Name
}
names(ALL.GENES) <- c("H. sapiens", "M. musculus", "R. norvegicus")
