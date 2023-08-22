#' @import data.table
#' @import Matrix
#' @importFrom stats fisher.test quantile xtabs



IBD <- function(ibd_data = "name.ibd", caco = "name.Rda", ...) {
  V1 <- V3 <- NULL
  ibd <- data.table::fread(ibd_data)
  ibd <- ibd[V1 %in% unlist(unique(caco)) & V3 %in% unlist(unique(caco))]
  class(ibd) <- append("IBD",class(ibd))
  return(ibd)
}



