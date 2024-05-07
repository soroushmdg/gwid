#' @import data.table
#' @import Matrix
#' @import shiny
#' @import plotly
#' @import lattice
#' @importFrom stats fisher.test quantile xtabs qbeta
#' @importFrom piggyback pb_download
#' @importFrom grid grid.polygon gpar





IBD <- function(ibd_data = "name.ibd", caco = "name.Rda", ...) {
  V1 <- V3 <- NULL
  ibd <- data.table::fread(ibd_data)
  ibd <- ibd[V1 %in% unlist(unique(caco)) & V3 %in% unlist(unique(caco))]
  class(ibd) <- append("IBD",class(ibd))
  return(ibd)
}

utils::globalVariables(c("caco_file_name", "data_folder_address"))



