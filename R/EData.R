EData <- methods::setClass("EData", methods::representation(
  taxon = "character", data = "list",
  experiment.meta = "data.table", gene.meta = "data.table",
  go = "data.table"
))
