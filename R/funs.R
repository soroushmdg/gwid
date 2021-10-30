#' Data Extraction of IBD and gds data
#'
#'Extract sample ID, IBD's for case-control subjects, snp ID, snp position, sample index, and smp.snp object to use in the pipeline.
#' @param gds_data CoreArray Genomic Data Structure (GDS) with hierarchical structure to store
#' multiple scalable structure array oriented data sets.
#' @param caco A List of case and control subjects.
#' @param ibd_data Data frame containing identity by descent segments in phased data. The data is output of refined_ibd
#' software package.
#' @param i chromosome number can be from 1 to 22
#' @param ... argument to pass to plots (if we decided to add graphs)
#'
#' @return data frames that are derived from inputs in order to overcome computational obstacle
#'
#' @examples new_gwid("chr1.gds","case_control.Rda","chr1.ibd", i = 1)
#' @export
new_gwas <- function(gds_data = "name.gds", caco = "name.Rda", ...) {
  substrR <- function(x, n) {
    substr(x, nchar(x) - n + 1, nchar(x))
  }
  if (substrR(caco, 4) != ".Rda") {
    csv2Rdat(caco)
    caco <- paste0(substr(caco, 1, nchar(caco) - 4), ".Rda")
  }
  load(caco)
  genoRA <- SNPRelate::snpgdsOpen(gds_data)
  smp.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "sample.id"))
  #IBD <- data.table::fread(ibd_data)
  #IBD <- IBD[V1 %in% unlist(unique(caco)) & V3 %in% unlist(unique(caco))]
  snp.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "snp.rs.id"))
  snp.pos <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "snp.position"))
  smp.indx <- which(smp.id %in% unique(unlist(caco)))
  smp.snp <- list()
  for (j in 1:length(caco)) {
    smp.snp[[j]] <- (gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "genotype"))[which(smp.id[smp.indx] %in% caco[[j]]), ])
    rownames(smp.snp[[j]]) <- smp.id[which(smp.id[smp.indx] %in% caco[[j]])]
    colnames(smp.snp[[j]]) <- snp.id
    smp.snp[[j]][smp.snp[[j]] == 3] <- NA
  }
  SNPRelate::snpgdsClose(genoRA)
  return(structure(list(smp.id = smp.id, snp.id = snp.id, snp.pos = snp.pos, smp.indx = smp.indx, smp.snp = smp.snp, caco = caco), class="snp_smp"))
}


#' @import data.table
#' @export
IBD <- function(ibd_data = "name.ibd" , caco = "name.Rda", ...){
  load(caco)
  ibd <- data.table::fread(ibd_data)
  ibd <- ibd[V1 %in% unlist(unique(caco)) & V3 %in% unlist(unique(caco))]
  class(ibd) <- append("IBD", class(ibd))
  return(ibd)
}


#' @export
print.snp_smp <- function(snp_smp,  ...) {
  print(paste("GWAS with", length(snp_smp$smp.id), "samples and", length(snp_smp$snp.id),"SNPs"))
  invisible(snp_smp)
}

#' @export
aggregate.snp_smp <- function(snp_smp, ...) {
  snp2 <- snps <- nas <- matrix(0L, nrow = length(snp_smp$snp.id), ncol = length(snp_smp$caco))
  for (j in 1:length(snp_smp$caco)) {
    nas[, j] <- as.integer(apply(is.na(snp_smp$smp.snp[[j]]), 2, sum))
    snps[, j] <- as.integer(apply(snp_smp$smp.snp[[j]], 2, sum, na.rm = T))
    snp2[, j] <- as.integer(apply(snp_smp$smp.snp[[j]] == 2, 2, sum, na.rm = T))
  }
  colnames(snp2) <- colnames(snps) <- colnames(nas) <- names(snp_smp$caco)
  # snp.chr <- rep(i, length(gwas$snp.id))
  output <- list(nas = nas, snp2 = snp2, snps = snps)
  class(output) <- append("gwas", class(output))
  return(output)
}


#' @export
plot.gwas <- function(statistics, data, type = "snps") {
  plot_general(statistics, data, type)
}

#' @export
phased <- function(x, ...){ UseMethod("phased") }

#' @export
phased.snp_smp <- function(snp_smp, phased_vcf, ...) {
  if (is.null(snp_smp) || is.null(phased_vcf)) stop("'snp_smp' and 'phased vcf' are needed")
  phased <- list()
  tmp <- data.table::fread(phased_vcf)
  tmp2 <- tmp[, colnames(tmp) %in% snp_smp$smp.id[snp_smp$smp.indx], with = F]
  phased[[1]] <- Matrix::Matrix(apply(apply(tmp2, 2, substr, 1, 1), 2, as.integer))
  phased[[2]] <- Matrix::Matrix(apply(apply(tmp2, 2, substr, 3, 3), 2, as.integer))
  names(phased) <- c("Hap.1", "Hap.2")
  class(phased) <- "phased"
  return(phased)
}



#' @export
new_gwid <- function(snp_smp, phased, IBD, ...) {
  mrk1 <- Mres <- LST <- ind <- list()
  INDX <- Matrix::Matrix(0, nrow = nrow(IBD), ncol = 6)
  colnames(INDX) <- c("Subj1", "Subj2", "start", "end", "hap1", "hap2")
  for (j in 1:length(snp_smp$caco)) {
    ind[[j]] <- which(IBD$V1 %in% snp_smp$caco[[j]] & IBD$V3 %in% snp_smp$caco[[j]])
    mrk1[[j]] <- LST[[j]] <- list()
    Mres[[j]] <- Matrix::Matrix(outer(IBD[ind[[j]], ]$V6, snp_smp$snp.pos, FUN = "<=") + outer(IBD[ind[[j]], ]$V7, snp_smp$snp.pos, FUN = ">=") - 1)
    for (k in 1:length(ind[[j]])) {
      INDX[ind[[j]][k], 1:2] <- which((snp_smp$smp.id[snp_smp$smp.indx]) %in% IBD[ind[[j]][k], c(1, 3)])
      INDX[ind[[j]][k], 3:4] <- which(snp_smp$snp.pos %in% IBD[ind[[j]][k], 6:7])
      INDX[ind[[j]][k], 5:6] <- unlist(IBD[ind[[j]][k], c(2, 4)])
      LST[[j]][[k]] <- Matrix::Matrix(phased[[INDX[ind[[j]][k], 5]]][INDX[ind[[j]][k], 3]:INDX[ind[[j]][k], 4], INDX[ind[[j]][k], 1]])
      mrk1[[j]][[k]] <- as.character(snp_smp$snp.id[INDX[ind[[j]][k], "start"]:INDX[ind[[j]][k], "end"]][which(LST[[j]][[k]][, 1] == 1)])
    }
  }
  names(ind) <- names(Mres) <- names(LST) <- names(snp_smp$caco)
  output <- list(mrk1 = mrk1, Mres = Mres, LST = LST, INDX = INDX, ind = ind)
  class(output) <- append("raw_gwid", class(output))
  return(output)
}

#' @export
aggregate.raw_gwid <- function(raw_gwid, snp_smp, IBD, ...) {
  res <- res1 <- IND <- Subj.id <- list()
  len <- NULL
  IND <- list()
  LST <- raw_gwid$LST
  Mres <- raw_gwid$Mres
  res <- matrix(0, nr = length(snp_smp$snp.pos), nc = length(snp_smp$caco))
  res1 <- matrix(0, nr = length(snp_smp$snp.pos), nc = length(snp_smp$caco))
  rownames(res1) <- snp_smp$snp.id
  for (j in 1:length(snp_smp$caco)) {
    IND[[j]] <- which(IBD$V1 %in% snp_smp$caco[[j]] & IBD$V3 %in% snp_smp$caco[[j]])
    res[, j] <- apply(Mres[[j]], 2, sum)
    temp <- table(unlist(raw_gwid$mrk1[[j]]))
    res1[names(temp), j] <- temp
  }
  names(IND) <- colnames(res) <- colnames(res1) <- names(snp_smp$caco)
  Subj.id <- (IBD[, c(1, 3)])
  output <- list(LST = LST, Mres = Mres, INDX = raw_gwid$INDX, Subj.id = Subj.id, IND = IND, res = res, res1 = res1)
  class(output) <- append("gwid", class(output))
  return(output)
}

usethis::use_pipe()



#' @export
plot.gwid <- function(statistics, data, type = "res") {
  plot_general(statistics, data, type)
}


#' @export
wind_base <- function(data, w) {
  data <- data[["Mres"]]
  if (w > ncol(data[[1]])) stop("window size must be smaller than colums size of the data")
  if (!(is.list(data))) stop("data must be of class list!")
  wind_sum <- matrix(0, nrow = (ncol(data[[1]]) - w + 1), ncol = length(data))
  for (j in 1:length(data)) {
    for (i in (1:(ncol(data[[1]]) - w + 1))) {
      wind_sum[i, j] <- sum(apply(data[[j]][, i:(w + i - 1), drop = F], 1, sum) == w)
    }
  }
  colnames(wind_sum) <- c("cases","case1","case2","cont1","cont2","cont3")
  output <- list(window = wind_sum)
  class(output) <- append("window", class(output))
  return(output)
}


#' @export
haplo_win_str <- function(raw_gwid, snp_smp, w = 10, snp_start, snp_end) {
  if (missing(snp_start)) {
    snp_start <- snp_smp$snp.pos[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_smp$snp.pos[length(snp_smp$snp.pos) - w + 1]
  }

  snp_indx <- (which(snp_smp$snp.pos >= snp_start & snp_smp$snp.pos <= snp_end))
  start_end_indx <- range(snp_indx)
  leni <- diff(start_end_indx) - w + 1
  if (w >= diff(start_end_indx)) {
    stop("window size should be smaller than number of snps")
  }
  lenj <- length(raw_gwid$Mres)
  structures <- vector(mode = "list", length = lenj)
  names(structures) <- names(raw_gwid$Mres)
  for (j in 1:lenj) {
    structures[[j]] <- vector(mode = "list", length = leni)
    for (i in 1:leni) {
      ibd_reg_count <- which(apply(raw_gwid$Mres[[j]][, snp_indx[i]:(w + snp_indx[i] - 1), drop = F], 1, sum) == w)
      for (k in 1:length(ibd_reg_count)) {
        structures[[j]][[i]][k] <- paste0(raw_gwid$LST[[j]][[ibd_reg_count[k]]][snp_indx[i]:(w + snp_indx[i] - 1)], collapse = "")
      }
    }
  }
  names(structures) <- c("cases","case1","case2","cont1","cont2","cont3")
  output <- list(structures = structures)
  return(output)
}

#' @export
haplo_win_frame <- function(raw_gwid, snp_smp, w = 10, snp_start, snp_end) {
  if (missing(snp_start)) {
    snp_start <- snp_smp$snp.pos[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_smp$snp.pos[length(snp_smp$snp.pos) - w + 1]
  }

  snp_indx <- (which(snp_smp$snp.pos >= snp_start & snp_smp$snp.pos <= snp_end))
  start_end_indx <- range(snp_indx)
  leni <- diff(start_end_indx) - w + 1
  if (w >= diff(start_end_indx)) {
    stop("window size should be smaller than number of snps")
  }
  lenj <- length(raw_gwid$Mres)
  haplo_win <- matrix(character(),ncol = 3)
  #names(structures) <- names(raw_gwid$Mres)
  for (j in 1:lenj) {
    #structures[[j]] <- vector(mode = "list", length = leni)
    for (i in 1:leni) {
      ibd_reg_count <- which(apply(raw_gwid$Mres[[j]][, snp_indx[i]:(w + snp_indx[i] - 1), drop = F], 1, sum) == w)
      for (k in 1:length(ibd_reg_count)) {
        structures<- paste0(raw_gwid$LST[[j]][[ibd_reg_count[k]]][snp_indx[i]:(w + snp_indx[i] - 1)], collapse = "")
        new <- c(caco = names(raw_gwid$Mres)[j], window = paste0("window",i),structure =  structures)
        haplo_win <- rbind(haplo_win,new)
      }
    }
  }
  #names(structures) <- c("cases","case1","case2","cont1","cont2","cont3")
  rownames(haplo_win) <- NULL
  output <- haplo_win
  return(output)
}

#' @export
haplo_win_freq <- function(hap_str) {
  freq <- matrix(unlist(lapply(lapply(unlist(hap_str[[1]],recursive = F),table),max)),ncol = length(hap_str[[1]]),byrow = F)
  colnames(freq) <- c("cases","case1","case2","cont1","cont2","cont3")
  return(freq)
}




#' @export
plot.window <- function(statistics, data, type = "window") {
  plot_general(statistics, data, type)
}


plot_general <- function(statistics, data, type = "snps") {
  if (type %in% c("snps", "nas", "snp2", "res","window")) {
    df <- tibble::as_tibble(cbind(snp.pos = data$snp.pos[1:length(statistics[[type]][,1])], statistics[[type]])) %>%
      tidyr::pivot_longer(!snp.pos, names_to = "case_control", values_to = "value")

    p <- df %>% ggplot2::ggplot(ggplot2::aes(x = snp.pos, y = value)) +
      ggplot2::geom_line(ggplot2::aes(color = case_control), size = .6) +
      ggplot2::scale_x_continuous("snp position", labels = scales::label_number_si()) +
      ggplot2::scale_y_continuous("sum of type in IBD regions") +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 10))

    fig <- plotly::ggplotly(p)
    fig
  }
}




csv2Rdat <- function(name = "", type = 1, rep = 3) {
  if (!file.exists(name)) message("File doesn't exist")
  set.seed(1)
  caco <- list()
  input.data <- read.csv(name)
  check.column <- function(v) {
    return(sum(c("1", "2") %in% as.integer(v)))
  }
  cacos <- names(which(apply(input.data, 2, check.column) == 2))

  for (i in 1:length(cacos)) {
    caco$cases <- input.data$SampleName[which(input.data[, cacos[i]] == 2)]
    cont.ind <- which(input.data[, cacos[i]] == 1)
    if (length(caco$cases) * rep > length(cont.ind)) caco$cases <- sample(caco$cases, trunc(length(cont.ind) / rep))
    conts <- sample(cont.ind, length(caco$cases) * rep)
    for (j in 1:rep) {
      caco[[paste0("cont", j)]] <- input.data$SampleName[conts[rep(1:rep, each = length(caco$cases)) == j]]
    }
    save(file = paste0(cacos[i], ".Rda"), caco)
  }
}

