# smaller data objects that are stored in package
# most of these objects used to be created on application runtime
# they are moved here to prevent computation during application load time and
# to simplify the loading process- ogan
print('package data')
library(magrittr)
devtools::load_all()
source('main/dependencies.R')
.DATA_PATH <- paste(DATADIR, 'DATA.HOLDER.light.rds', sep='/')


# ontology data ------
# backed up nathaniels files at cosmos with the rest of the data. keeping these
# to remember the original paths if needed
# narrator: they were needed
# ONTOLOGIES <- data.table::fread("/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED.TSV")
# ONTOLOGIES.DEFS <- data.table::fread("/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED_DEF.TSV")


fixOntoGenes <- function() {
  lapply(getConfig(key = "taxa")$core, function(taxa) {
    matches <- do.call(
      rbind,
      stringr::str_match_all(
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


#ONTOLOGIES <- data.table::fread(file.path(FREEZEDIR,"CronGemmaDump/Ontology/Ontology_Dump_MERGED.TSV"))
#ONTOLOGIES.DEFS <- data.table::fread(file.path(FREEZEDIR,"CronGemmaDump/Ontology/Ontology_Dump_MERGED_DEF.TSV"))

# re-using nathaniel's files again as they are updated monthly
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


if (file.exists(.DATA_PATH)) {
  DATA.HOLDER <- readRDS(.DATA_PATH)
} else {
  stop("Couldn't find DATA.HOLDER, run generate/fbm.R first.")
}
  
DATA.HOLDER$artificial <- NULL
  
ONTOLOGIES.DEFS <- fixOntoGenes()



# simplified ontology defs removes the ontology scope, resolves ambiguous terms

# identify duplicates and disambiguate based on majority
dups <- ONTOLOGIES.DEFS %>% 
  dplyr::select(Node_Long,Definition) %>% 
  unique %$% Node_Long %>% {.[duplicated(.)]}

# a bunch of heuristics to disambiguate
dups_resolved <- dups %>% parallel::mclapply(function(x){
  defs = ONTOLOGIES.DEFS %>% 
    dplyr::filter(Node_Long == x) %$% 
    Definition %>% as.character() %>% 
    tolower() %>% 
    table %>% 
    sort(decreasing = TRUE)
  if(length(defs)== 1 || defs[1] == defs[2]){
    # if a majority is not present try to extract the context and listen to the original source
    context <- x %>% stringr::str_extract("(?<=[a-z]/)[A-Z]*?(?=_)")
    if(!is.na(context) && context == 'DOID'){
      context = 'DO'
    }
    disambig <- ONTOLOGIES.DEFS %>% dplyr::filter(Node_Long == x & OntologyScope == context) %$% Definition 
    if(length(disambig) == 1){
      return(as.character(disambig))
    } else{
      # finally, if all fails, go with the longest non-number definition 
      possible =  ONTOLOGIES.DEFS %>% dplyr::filter(Node_Long == x) %$% Definition
      char_count <- possible %>% stringr::str_replace_all("[0-9]",'') %>% nchar()
      return(as.character(possible[which.max(char_count)]))
    }
  } else{
    # if a majority is present, go with it
    return(names(defs[1]))
  }
},mc.cores = 16) %>% unlist()
names(dups_resolved) = dups


SIMPLIFIED.ONTOLOGY.DEFS <- ONTOLOGIES.DEFS %>% dplyr::select(Node_Long,Definition) %>% unique()

SIMPLIFIED.ONTOLOGY.DEFS %<>% dplyr::filter(!Node_Long %in% names(dups_resolved))
SIMPLIFIED.ONTOLOGY.DEFS %<>% rbind(
  data.table(Node_Long = dups_resolved %>% names,
             Definition = dups_resolved)
)



usethis::use_data(ONTOLOGIES,overwrite = TRUE)
usethis::use_data(ONTOLOGIES.DEFS,overwrite = TRUE)
usethis::use_data(SIMPLIFIED.ONTOLOGY.DEFS,overwrite = TRUE)


# CACHE.BACKGROUND ------------
CACHE.BACKGROUND <- lapply(names(DATA.HOLDER), function(x){
  ONTOLOGIES = ONTOLOGIES
  # ONTOLOGIES.DEFS = ONTOLOGIES.DEFS
  precomputeTags(x)
}) %>%
  setNames(names(DATA.HOLDER))

usethis::use_data(CACHE.BACKGROUND,overwrite = TRUE)

# NULLS --------------

NULLS <- generate_nulls('score_contrast')

usethis::use_data(NULLS,overwrite = TRUE)


# ALL.GENES -------
ALL.GENES <- list()
for (taxon in names(DATA.HOLDER)) {
  ALL.GENES[[taxon]] <- DATA.HOLDER[[taxon]]@gene.meta$gene.Name
}
names(ALL.GENES) <- c("H. sapiens", "M. musculus", "R. norvegicus")

usethis::use_data(ALL.GENES,overwrite = TRUE)


# taxonomy information -------
tax_data = homologene::taxData[match(c(9606, 10090, 10116),homologene::taxData$tax_id),]

TAX.DATA= data.frame(id = tax_data$tax_id,
                     common_names = c("human", "mouse", "rat"),
                     names = tax_data$name_txt
                     )

usethis::use_data(TAX.DATA,overwrite = TRUE)


# filters ---------------

child_filters = jsonlite::fromJSON(readLines('data-raw/filters.json'))
text_filters = jsonlite::fromJSON(readLines('data-raw/text_filters.json'))

filters = list(base_filter = c(get_parents(child_filters$base_filter),text_filters$base_filter),
               val_filter = c(get_parents(child_filters$val_filter),text_filters$val_filter),
               universal_filter = c(get_parents(child_filters$universal_filter),text_filters$universal_filter))


usethis::use_data(filters,overwrite = TRUE)
