#' runs of homozygosity
#'
#' @param phase object of phase
#'
#' @param ... other variables
#'
#' @return runs of homozygosity data table or matrix
#' @export
roh <- function(phase, ...) {
  UseMethod("roh")
}

#' runs of homozygosity
#'
#' @param phase An object of class phase. Output of \code{build_phase} function
#'
#' @param gwas object of class gwas
#' @param w window size
#' @param fun an aggregate function. either  \dQuote{sum} or  \dQuote{mean}
#' @param snp_start select starting position of snps.
#' @param snp_end select ending position of snps.
#' @param roh_mat return roh as matrix
#' @param ... other variables
#'
#' @return the output will be a result_snps (data.table) object including 3 columns
#' including, \dQuote{snp_pos}, \dQuote{case_control}, and \dQuote{value}
#'
#'
#' @export
roh.phase <- function(phase, gwas, w = 10, fun = c("sum", "mean"), snp_start, snp_end, roh_mat = FALSE, ...) {
  if (missing(phase) | missing(gwas)) {
    stop("Please provide function arguments")
  }
  fun <- match.arg(fun)
  if (missing(snp_start)) {
    snp_start <- gwas$snp.pos[1]
  }
  if (missing(snp_end)) {
    snp_end <- gwas$snp.pos[length(gwas$snp.pos)]
  }
  snp_indx <- which(gwas$snp.pos >= snp_start & gwas$snp.pos <= snp_end)
  myphase <- lapply(phase, function(x) x[snp_indx, ])
  haps_diff <- myphase$Hap.1 == myphase$Hap.2
  haps_diff <- haps_diff * 1
  ROH <- haps_diff[1:(nrow(haps_diff) - (w - 1)), ]
  for (i in 1:(w - 1)) {
    ROH <- ROH + haps_diff[(i + 1):(nrow(haps_diff) - (w - 1) + i), ]
  }
  if (roh_mat == TRUE) {
    return(ROH)
  }
  snp_pos <- gwas$snp.pos[snp_indx][1:nrow(ROH)]
  Roh_stat <- matrix(0, nrow(ROH), length(gwas[[which(unlist(lapply(gwas, inherits, "caco")))]]))
  for (i in 1:ncol(Roh_stat)) {
    temp <- gwas[[which(unlist(lapply(gwas, inherits, "caco")))]][[i]]
    roh1 <- ROH[, temp]
    Roh_stat[, i] <- apply(roh1 == w, 1, fun)
  }
  colnames(Roh_stat) <- names(gwas[[which(unlist(lapply(gwas, inherits, "caco")))]])
  Roh_stat <- data.table::as.data.table(cbind(Roh_stat, snp_pos = snp_pos))
  roh_stat <- data.table::melt(Roh_stat, id.vars = "snp_pos", variable.name = "case_control", value.name = "value")
  class(roh_stat) <- append("result_snps", class(roh_stat))
  return(roh_stat)
}



#' mcnemar test
#'
#' @param roh roh as class result_snp
#'
#' @param ... other variables
#'
#' @export
mcnemar_test <- function(roh, ...) {
  UseMethod("mcnemar_test")
}


#' mcnemar test
#'
#' @param roh An object of class result_snps (output of function roh with fun=sum)
#'
#' @param reference reference group of subjects in which we want to perform fisher test.
#' @param w window size
#' @param ... other variables
#'
#' @export
mcnemar_test.result_snps <- function(roh = "object of class result_snps (output of function roh with fun=sum)", reference, w = 10, ...) {
  if (missing(reference)) {
    stop("Please provide a reference e.x. 'case1' ")
  }
  snp_pos <- case_control <- NULL
  Roh_stat <- data.table::dcast(roh, snp_pos ~ case_control, value.var = "value")
  snp_pos <- Roh_stat$snp_pos
  Roh_stat <- as.matrix(Roh_stat[, !"snp_pos", with = FALSE])

  # ref_num: identify reference
  ref_num <- which(colnames(Roh_stat) == reference)
  McNemar <- matrix(0, nrow(Roh_stat), ncol(Roh_stat))
  McNemar_smooth <- matrix(0, nrow(Roh_stat) - w + 1, ncol(Roh_stat))
  for (i in 1:ncol(Roh_stat)) {
    McNemar[, i] <- (Roh_stat[, ref_num] - Roh_stat[, i])^2 / (Roh_stat[, ref_num] + Roh_stat[, i])
    McNemar_smooth[, i] <- RcppRoll::roll_mean(x = McNemar[, i], n = w, align = "right")
  }
  colnames(McNemar_smooth) <- colnames(Roh_stat)
  McNemar_smooth <- data.table::as.data.table(cbind(McNemar_smooth, snp_pos = snp_pos[1:nrow(McNemar_smooth)]))
  McNemar <- data.table::melt(McNemar_smooth, id.vars = "snp_pos", variable.name = "case_control", value.name = "value")
  output <- McNemar
  class(output) <- append("result_snps", class(output))
  return(output)
}


#' mcnemar permutation
#'
#' @param mcnemar macnemar test output
#'
#' @param ... other variables
#'
#' @export
mcnemar_test_permut <- function(mcnemar, ...) {
  UseMethod("mcnemar_test_permut")
}


#' mcnemar permutation test
#'
#' @param mcnemar macnemar test output
#'
#' @param roh_mat roh matrix
#' @param gwas gwas
#' @param nperm number of permutation
#' @param reference reference group
#' @param w window
#' @param ... other variables
#'
#' @export
mcnemar_test_permut.result_snps <- function(mcnemar = "object of class result_snps (output of function mcnemar_test with fun=sum)",
                                            roh_mat = "output of roh function when roh_mat = TRUE",
                                            gwas = "object of class gwas",
                                            nperm = 1000,
                                            reference = "cases",
                                            w, ...) {
  snp_pos <- case_control <- NULL
  mymcnemar1 <- data.table::dcast(mcnemar, snp_pos ~ case_control, value.var = "value") # [,snp_pos:=NULL])
  snp_pos <- mymcnemar1$snp_pos
  mymcnemar1 <- as.matrix(mymcnemar1[, snp_pos := NULL])
  mc_perm <- vector(mode = "list", length = nperm)
  for (i in 1:length(mc_perm)) {
    roh_perm <- roh_mat
    colnames(roh_perm) <- sample(colnames(roh_mat))
    Roh_stat <- matrix(0, nrow = nrow(roh_mat), ncol = ncol(mymcnemar1))

    temp_ind <- lapply(gwas[["caco"]], function(x) {
      which(colnames(roh_perm) %in% x)
    })

    Roh_stat <- sapply(temp_ind, function(x) {
      roh1 <- roh_perm[, x]
      rowSums(roh1 == w)
    })
    ref_num <- which(colnames(Roh_stat) == reference)
    McNemar_smooth <- sapply((1:ncol(Roh_stat)), function(x) {
      McNemar <- (Roh_stat[, ref_num] - Roh_stat[, x])^2 / (Roh_stat[, ref_num] + Roh_stat[, x])
      RcppRoll::roll_mean(x = McNemar, n = w, align = "right")
    })
    colnames(McNemar_smooth) <- colnames(Roh_stat)
    mc_perm[[i]] <- McNemar_smooth
  }
  mc_perm <- do.call(what = "cbind", mc_perm)
  pval <- matrix(0, nrow = nrow(mymcnemar1), ncol = ncol(mymcnemar1))
  for (i in colnames(mymcnemar1)) {
    pval[, which(colnames(mymcnemar1) %in% i)] <- apply(apply(mc_perm[, which(colnames(mc_perm) == i)], 2, ">", mymcnemar1[, i]), 1, mean) #
  }
  colnames(pval) <- colnames(mymcnemar1)
  pval[, which(colnames(pval) == reference)] <- 1
  pval <- data.table::as.data.table(pval)
  pval[, snp_pos := snp_pos]
  pval <- data.table::melt(pval, id.vars = "snp_pos", variable.name = "case_control", value.name = "value")
  class(pval) <- append("test_snps", class(pval))
  return(pval)
}
