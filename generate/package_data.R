# smaller data objects that are stored in package
# most of these objects used to be created on application runtime
# they are moved here to prevent computation during application load time and
# to simplify the loading process- ogan

devtools::load_all()
.DATA_PATH <- paste(DATADIR, 'DATA.HOLDER.light.rds', sep='/')


# ontology data ------
# backed up nathaniels files at cosmos with the rest of the data. keeping these
# to remember the original paths if needed
# ONTOLOGIES <- data.table::fread("/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED.TSV")
# ONTOLOGIES.DEFS <- data.table::fread("/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED_DEF.TSV")

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


ONTOLOGIES <- data.table::fread(file.path(FREEZEDIR,"CronGemmaDump/Ontology/Ontology_Dump_MERGED.TSV"))
ONTOLOGIES.DEFS <- data.table::fread(file.path(FREEZEDIR,"CronGemmaDump/Ontology/Ontology_Dump_MERGED_DEF.TSV"))


ONTOLOGIES[, c("ChildNode", "ParentNode")] <- NULL
ONTOLOGIES.DEFS$Node <- NULL

ONTOLOGIES$ChildNode_Long <- ONTOLOGIES$ChildNode_Long %>% as.factor()
ONTOLOGIES$ParentNode_Long <- ONTOLOGIES$ParentNode_Long %>% as.factor()
ONTOLOGIES$RelationType <- ONTOLOGIES$RelationType %>% as.factor()
ONTOLOGIES$OntologyScope <- ONTOLOGIES$OntologyScope %>% as.factor()

ONTOLOGIES.DEFS$Node_Long <- ONTOLOGIES.DEFS$Node_Long %>% as.factor()
ONTOLOGIES.DEFS$OntologyScope <- ONTOLOGIES.DEFS$OntologyScope %>% as.factor()

if (file.exists(.DATA_PATH)) {
  DATA.HOLDER <- readRDS(.DATA_PATH)
} else {
  stop("Couldn't find DATA.HOLDER, run generate/fbm.R first.")
}
  
DATA.HOLDER$artificial <- NULL
  
ONTOLOGIES.DEFS <- fixOntoGenes()


usethis::use_data(ONTOLOGIES,overwrite = TRUE)
usethis::use_data(ONTOLOGIES.DEFS,overwrite = TRUE)


# CACHE.BACKGROUND ------------
CACHE.BACKGROUND <- lapply(names(DATA.HOLDER), precomputeTags) %>%
  setNames(names(DATA.HOLDER))

usethis::use_data(CACHE.BACKGROUND,overwrite = TRUE)

# NULLS --------------

NULLS <- lapply(names(DATA.HOLDER), function(taxa) {
  tryCatch(readRDS(paste0(paste(DATADIR, 'updated_nulls/', sep='/'), taxa, ".rds")) %>%
             .[, .(rn, score.mean, score.sd)], error = function(e) NULL)
}) %>%
  `names<-`(names(DATA.HOLDER)) %>%
  {
    Filter(Negate(is.null), .)
  }

usethis::use_data(NULLS,overwrite = TRUE)


# ALL.GENES -------
ALL.GENES <- list()
for (taxon in names(DATA.HOLDER)) {
  ALL.GENES[[taxon]] <- DATA.HOLDER[[taxon]]@gene.meta$gene.Name
}
names(ALL.GENES) <- c("H. sapiens", "M. musculus", "R. norvegicus")

usethis::use_data(ALL.GENES,overwrite = TRUE)


# taxonomy information -------
tax_data = homologene::taxData %>% filter(tax_id %in%  c(9606, 10090, 10116))

TAX.DATA= data.frame(id = tax_data$tax_id,
                     )

