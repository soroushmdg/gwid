
#' Genome Wide Identity by Descent (GWID)
#'
#' @param gds_data CoreArray Genomic Data Structure (GDS) with hierarchical structure to store
#' multiple scalable structure array oriented data sets.
#'
#' @param caco A List of case and control subjects.
#' @param ibd_data Data frame containing identity by descent segments in phased data. The data is output of refined_ibd
#' software package.
#' @param vcf_data Phased genotype data in vcf format.
#' @param i chromosome number can be from 1 to 22
#' @param ... argument to pass to plots (if we decided to add graphs)
#'
#' @return data frames that are derived from inputs in order to overcome computational obstacle.
#' @export
#'
#' @examples
gwid <- function(gds_data, caco, ibd_data, vcf_data ,i,...) {
  V1 <- V2 <- V3 <- V4 <- V5 <- V6 <- V7 <- V8<- V9 <- NULL
IBD <- genoRA <- snp.id <- snp.pos <- smp.id <- smp.indx <- smp.snp <- list()
  genoRA <- SNPRelate::snpgdsOpen(gds_data)
  smp.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "sample.id"))
  IBD <- data.table::fread(ibd_data)[V1 %in% unlist(unique(caco)) & V3 %in% unlist(unique(caco))]
  snp.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "snp.rs.id"))
  snp.pos <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "snp.position"))
  smp.indx <- which(smp.id %in% unique(unlist(IBD[, c(1, 3)])))
  smp.snp <- list()
  for (j in 1:length(caco)) {
    smp.snp[[j]] <- (gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "genotype"))[which(smp.id[smp.indx] %in% caco), ])
    rownames(smp.snp[[j]]) <- smp.id[which(smp.id[smp.indx] %in% caco[[j]])]
    colnames(smp.snp[[j]]) <- snp.id
    smp.snp[[j]][smp.snp[[j]] == 3] <- NA
  }
  SNPRelate::snpgdsClose(genoRA)


  snp2 <- snps <- nas <- matrix(0L, nrow = length(snp.id), ncol = length(caco))
  for (j in 1:length(caco)) {
    nas[, j] <- as.integer(apply(is.na(smp.snp[[j]]), 2, sum))
    snps[, j] <- as.integer(apply(smp.snp[[j]], 2, sum, na.rm = T))
    snp2[, j] <- as.integer(apply(smp.snp[[j]] == 2, 2, sum, na.rm = T))
  }
  colnames(snp2) <- colnames(snps) <- colnames(nas) <- names(caco)

snp.chr <- rep(i, sapply(snp.id, length))
#snp.id <- unlist(snp.id)
#snp.pos <- unlist(snp.pos)


ibd <- IBD
  phased <- list()
  IBD <- ibd
  tmp <- data.table::fread(vcf_data)
  tmp2 <- tmp[, colnames(tmp) %in% smp.id[smp.indx], with = F]
  phased[[1]] <- Matrix::Matrix(apply(apply(tmp2, 2, substr, 1, 1), 2, as.integer))
  phased[[2]] <- Matrix::Matrix(apply(apply(tmp2, 2, substr, 3, 3), 2, as.integer))
  names(phased) <- c("Hap.1", "Hap.2")

rm(list = "ibd")

#load(paste0("../../snp.smp-RA.withmap.Rda"))
#

  #load(paste0("../../phased-IBD-RA.withmap-chr", i, ".Rda")) # doesn't run
  mrk1 <- Mres <- LST <- ind <- list()
  INDX <- Matrix::Matrix(0, nrow = nrow(IBD), ncol = 6)
  colnames(INDX) <- c("Subj1", "Subj2", "start", "end", "hap1", "hap2")
  for (j in 1:length(caco)) {
    ind[[j]] <- which(IBD$V1 %in% caco[[j]] & IBD$V3 %in% caco[[j]])
    mrk1[[j]] <- LST[[j]] <- list()
    Mres[[j]] <- Matrix::Matrix(outer(IBD[ind[[j]], ]$V6, snp.pos, FUN = "<=") + outer(IBD[ind[[j]], ]$V7, snp.pos, FUN = ">=") - 1)
    for (k in 1:length(ind[[j]])) {
      INDX[ind[[j]][k], 1:2] <- which((smp.id[smp.indx]) %in% IBD[ind[[j]][k], c(1, 3)])
      INDX[ind[[j]][k], 3:4] <- which(snp.pos %in% IBD[ind[[j]][k], 6:7])
      INDX[ind[[j]][k], 5:6] <- unlist(IBD[ind[[j]][k], c(2, 4)])
      LST[[j]][[k]] <- Matrix::Matrix(phased[[INDX[ind[[j]][k], 5]]][INDX[ind[[j]][k], 3]:INDX[ind[[j]][k], 4], INDX[ind[[j]][k], 1]])
      mrk1[[j]][[k]] <- as.character(snp.id[INDX[ind[[j]][k], "start"]:INDX[ind[[j]][k], "end"]][which(LST[[j]][[k]][, 1] == 1)])
    }
  }
  names(Mres) <- names(LST) <- names(caco)


len <- NULL
  IND <- list()
  res <- matrix(0, nrow = length(snp.pos), ncol = length(caco))
  res1 <- matrix(0, nrow = length(snp.pos), ncol = length(caco))
  rownames(res1) <- snp.id
  for (j in 1:length(caco)) {
    IND[[j]] <- which(IBD$V1 %in% caco[[j]] & IBD$V3 %in% caco[[j]])
    res[, j] <- apply(Mres[[j]], 2, sum)
    temp <- table(unlist(mrk1[[j]]))
    res1[names(temp), j] <- temp
  }
  names(IND) <- colnames(res) <- colnames(res1) <- names(caco)
  rownames(res1) <- NULL
  ind <- which(IND[[1]] %in% unlist(IND[2:3]))
  LST[[1]] <- LST[[1]][-ind]
  Mres[[1]] <- Mres[[1]][-ind, ]
  IND[[1]] <- IND[[1]][-ind]
  Subj.id <- as.matrix(IBD[, c(1, 3)])
  len <- c(len, nrow(IBD))


Subj.id <- as.data.frame(do.call(rbind, Subj.id))
len <- c(0, cumsum(len)[-length(len)])
if (i!=0){
for (j in 1:length(IND)) IND[[j]] <- IND[[j]] + len}



indx <- list()
shift.adj <- which((snp.chr[-length(snp.chr)] - snp.chr[-1]) != 0)
shift.adj <- c(0, shift.adj[-22])

  indx <- cbind(INDX[, 1:2], "chr" = 0, INDX[, 3:6])
  indx[indx[, 1] != 0, 3] <- i


ind <- indx
ind[] <- 0
  INDX <- ind
  INDX <- indx
  INDX <- do.call(rbind, INDX)
  INDX[INDX[, 3] == i, 4:5] <- INDX[INDX[, 3] == i, 4:5] + shift.adj[i]
  Mres1 <- phased
  names(Mres1[[i]]) <- c("Hap1", "Hap2")
  return(snp.pos = snp.pos, snp.id = snp.id, smp.id = smp.id, smp.indx = smp.indx,
         snp.chr = snp.chr, nas = nas, snps = snps, snp2 = snp2, phased = phased,
         LST = LST, Mres = Mres, INDX = INDX, IND = IND, Subj.id = Subj.id, res = res, res1 = res1, Mres1 = Mres1)

}
