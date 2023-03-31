devtools::load_all()

mGraph <- igraph::simplify(igraph::graph_from_data_frame(ONTOLOGIES[, .(as.character(ChildNode_Long), as.character(ParentNode_Long))]))

# blacklist definition

get_parents = function(terms){
  terms %>% lapply(function(x){
    igraph::subcomponent(mGraph, x, "out") %>% names
  }) %>% unlist %>% unique
}


universal_filter = 
  c(
    "http://purl.obolibrary.org/obo/BFO_0000023",
    "http://purl.obolibrary.org/obo/CHEBI_23888"
  )

val_filter = 
  c("http://purl.obolibrary.org/obo/OBI_0100026",
    "http://purl.obolibrary.org/obo/CHEBI_60004",
    "http://purl.obolibrary.org/obo/BFO_0000023",
    "http://purl.obolibrary.org/obo/DOID_4",
    "http://www.ebi.ac.uk/efo/EFO_0000408",
    "http://purl.obolibrary.org/obo/BFO_0000030",
    "http://purl.obolibrary.org/obo/CHEBI_23367",
    "http://purl.obolibrary.org/obo/CL_0000000",
    "http://purl.obolibrary.org/obo/UBERON_0000062",
    "http://purl.obolibrary.org/obo/CHEBI_25367",
    "http://purl.obolibrary.org/obo/CHEBI_24432",
    "http://purl.obolibrary.org/obo/CHEBI_51086",
    "http://www.ebi.ac.uk/efo/EFO_0000324",
    "http://purl.obolibrary.org/obo/CHEBI_23888",
    "http://www.ebi.ac.uk/efo/EFO_0001899",
    "http://purl.obolibrary.org/obo/OBI_0000025",
    "http://purl.obolibrary.org/obo/CHEBI_25212",
    "http://purl.obolibrary.org/obo/UBERON_0000064",
    "http://www.ebi.ac.uk/efo/EFO_0005168",
    "http://purl.obolibrary.org/obo/UBERON_0000477",
    "http://purl.obolibrary.org/obo/UBERON_0002530",
    "http://purl.obolibrary.org/obo/OBI_0000047",
    "http://purl.obolibrary.org/obo/UBERON_0000479",
    "http://www.ebi.ac.uk/efo/EFO_0001824",
    "http://www.ebi.ac.uk/efo/EFO_0001461",
    "http://www.ebi.ac.uk/efo/EFO_0002694",
    "http://www.ebi.ac.uk/efo/EFO_0004425",
    "http://purl.obolibrary.org/obo/OBI_0000220",
    "http://purl.obolibrary.org/obo/OBI_0000025",
    "http://purl.obolibrary.org/obo/RO_0002577",
    "http://gemma.msl.ubc.ca/ont/TGEMO_00003",
    "http://purl.obolibrary.org/obo/BFO_0000019",
    "http://purl.obolibrary.org/obo/PATO_0000049",
    "http://purl.obolibrary.org/obo/PATO_0000396",
    "http://purl.obolibrary.org/obo/PATO_0000395",
    "http://purl.obolibrary.org/obo/PATO_0000001",
    "http://purl.obolibrary.org/obo/UBERON_0000063",
    "http://purl.obolibrary.org/obo/UBERON_0004923",
    "http://purl.obolibrary.org/obo/DOID_225",
    "http://purl.obolibrary.org/obo/DOID_0060035",
    "http://purl.obolibrary.org/obo/UBERON_0000475"
  )

base_filter = c("http://gemma.msl.ubc.ca/ont/TGEMO_00000")



child_filters = list(base_filter = base_filter,
                     val_filter = val_filter,
                     universal_filter = universal_filter)

jsonlite::toJSON(child_filters,pretty = TRUE) %>% cat(file = 'data-raw/child_filters.json')

base_text_filter= c()
val_text_filter = c()
universal_text_filter = c("http://purl.obolibrary.org/obo/DOID_556")



text_filters = list(base_filter = base_text_filter,
                    val_filter = val_text_filter,
                    universal_filter = universal_text_filter)


jsonlite::toJSON(child_filters,pretty = TRUE) %>% cat(file = 'data-raw/text_filters.json')
