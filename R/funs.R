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
new_gwas <- function(gds_data = "name.gds", caco = "name.Rda", ibd_data = "name.ibd", ...) {
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
  IBD <- data.table::fread(ibd_data)
  IBD <- IBD[V1 %in% unlist(unique(caco)) & V3 %in% unlist(unique(caco))]
  snp.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "snp.rs.id"))
  snp.pos <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "snp.position"))
  smp.indx <- which(smp.id %in% unique(unlist(IBD[, c(1, 3)])))
  smp.snp <- list()
  for (j in 1:length(caco)) {
    smp.snp[[j]] <- (gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "genotype"))[which(smp.id[smp.indx] %in% caco[[j]]), ])
    rownames(smp.snp[[j]]) <- smp.id[which(smp.id[smp.indx] %in% caco[[j]])]
    colnames(smp.snp[[j]]) <- snp.id
    smp.snp[[j]][smp.snp[[j]] == 3] <- NA
  }
  SNPRelate::snpgdsClose(genoRA)
  return(structure(list(smp.id = smp.id, IBD = IBD, snp.id = snp.id, snp.pos = snp.pos, smp.indx = smp.indx, smp.snp = smp.snp, caco = caco), class="gwas"))
}

#' @import data.table
#' @export
print.gwas <- function(gwas, ...) {
  print(paste("GWAS with", length(gwas$smp.id), "samples and", nrow(gwas$IBD),"IBD regions over", length(gwas$snp.id),"SNPs"))
  invisible(gwas)
}

#' @export
aggregate.gwas <- function(gwas, ...) {
  snp2 <- snps <- nas <- matrix(0L, nrow = length(gwas$snp.id), ncol = length(gwas$caco))
  for (j in 1:length(gwas$caco)) {
    nas[, j] <- as.integer(apply(is.na(gwas$smp.snp[[j]]), 2, sum))
    snps[, j] <- as.integer(apply(gwas$smp.snp[[j]], 2, sum, na.rm = T))
    snp2[, j] <- as.integer(apply(gwas$smp.snp[[j]] == 2, 2, sum, na.rm = T))
  }
  colnames(snp2) <- colnames(snps) <- colnames(nas) <- names(gwas$caco)
  # snp.chr <- rep(i, length(gwas$snp.id))
  return(list(nas = nas, snp2 = snp2, snps = snps))
}

#' @export
phased <- function(x, ...){ UseMethod("phased") }

#' @export
phased.gwas <- function(gwas, phased_vcf, ...) {
  if (is.null(gwas) || is.null(phased_vcf)) stop("'gwas' and 'phased vcf' are needed")
  phased <- list()
  tmp <- data.table::fread(phased_vcf)
  tmp2 <- tmp[, colnames(tmp) %in% gwas$smp.id[gwas$smp.indx], with = F]
  phased[[1]] <- Matrix::Matrix(apply(apply(tmp2, 2, substr, 1, 1), 2, as.integer))
  phased[[2]] <- Matrix::Matrix(apply(apply(tmp2, 2, substr, 3, 3), 2, as.integer))
  names(phased) <- c("Hap.1", "Hap.2"); class(phased)="phased"
  return(phased)
}

# #' @export
# gwid <- structure(list(), class = "gwid")

# #' @export
# new_gwid <- function(x, ...){UseMethod("new_gwid") }

#' @export
new_gwid <- function(gwas = "", phased= "", ...){
  mrk1 <- Mres <- LST <- ind <- list()
  INDX <- Matrix::Matrix(0, nrow = nrow(gwas$IBD), ncol = 6)
  colnames(INDX) <- c("Subj1", "Subj2", "start", "end", "hap1", "hap2")
  for (j in 1:length(gwas$caco)) {
    ind[[j]] <- which(gwas$IBD$V1 %in% gwas$caco[[j]] & gwas$IBD$V3 %in% gwas$caco[[j]])
    mrk1[[j]] <- LST[[j]] <- list()
    Mres[[j]] <- Matrix::Matrix(outer(gwas$IBD[ind[[j]], ]$V6, gwas$snp.pos, FUN = "<=") + outer(gwas$IBD[ind[[j]], ]$V7, gwas$snp.pos, FUN = ">=") - 1)
    for (k in 1:length(ind[[j]])) {
      INDX[ind[[j]][k], 1:2] <- which((gwas$smp.id[gwas$smp.indx]) %in% gwas$IBD[ind[[j]][k], c(1, 3)])
      INDX[ind[[j]][k], 3:4] <- which(gwas$snp.pos %in% gwas$IBD[ind[[j]][k], 6:7])
      INDX[ind[[j]][k], 5:6] <- unlist(gwas$IBD[ind[[j]][k], c(2, 4)])
      LST[[j]][[k]] <- Matrix::Matrix(phased[[INDX[ind[[j]][k], 5]]][INDX[ind[[j]][k], 3]:INDX[ind[[j]][k], 4], INDX[ind[[j]][k], 1]])
      mrk1[[j]][[k]] <- as.character(gwas$snp.id[INDX[ind[[j]][k], "start"]:INDX[ind[[j]][k], "end"]][which(LST[[j]][[k]][, 1] == 1)])
    }
  }
  names(Mres) <- names(LST) <- names(gwas$caco)
  return(structure(list(mrk1 = mrk1, Mres = Mres, LST = LST, INDX = INDX, ind = ind),class = "gwid"))
}

#' @export
aggregate.gwid <- function(gwid, new_gwas , ...) {
   res <- res1 <- IND <- Subj.id <- list();
  len <- NULL
   IND <- list()
     LST <- gwid$LST
     Mres <- gwid$Mres
     res <- matrix(0, nr=length(new_gwas$snp.pos), nc=length(new_gwas$caco));
     res1<- matrix(0, nr=length(new_gwas$snp.pos), nc=length(new_gwas$caco));
     rownames(res1) <- new_gwas$snp.id
     for (j in 1:length(new_gwas$caco)) {
       IND[[j]] <- which(new_gwas$IBD$V1%in%new_gwas$caco[[j]] & new_gwas$IBD$V3%in%new_gwas$caco[[j]])
       res[,j] <- apply(Mres[[j]],2,sum)
       temp <- table(unlist(gwid$mrk1[[j]]));
       res1[names(temp),j] <- temp
     }
     names(IND) <- colnames(res) <- colnames(res1) <- names(new_gwas$caco);
     rownames(res1) <- NULL
     ind <- which(IND[[1]]%in%unlist(IND[2:3]));
     LST[[1]] <- LST[[1]][-ind];
     Mres[[1]] <- Mres[[1]][-ind,];
     IND[[1]] <- IND[[1]][- ind];
     Subj.id <- (new_gwas$IBD[,c(1,3)]);
     len <- c(len,nrow(new_gwas$IBD))
     #save("LST","Mres","INDX",file=paste0("../../RA.withmap/chr",i,".Rda"), version = 2)
     #return(LST = LST, Mres = Mres, INDX = INDX, subj.id = subj.id)

   #load("../../genotype-RA.withmap.Rda");
   #Subj.id <- as.data.frame(subj.id);
   len <- c(0,cumsum(len)[-length(len)]) # len is going to be zero when working with one file (one chromosome)
    #for (j in 1:length(IND)) IND[[j]] <- IND[[j]] + len[i]
   #save("IND","Subj.id","res","res1","snps","nas","snp2","snp.chr","snp.id","snp.pos", file="../../RA.withmap/IBD-SNP-RA.withmap.Rda", version = 2)
   return(list(LST = LST, Mres = Mres, INDX = gwid$INDX, Subj.id = Subj.id, IND= IND, res= res, res1 = res1))
}




csv2Rdat <- function(name="", type=1, rep=3) {
  if (!file.exists(name)) message("File doesn't exist")
  set.seed(1); caco <- list()
  input.data <- read.csv(name)
  check.column <- function(v) {return(sum(c("1","2") %in% as.integer(v)))}
  cacos <- names(which(apply(input.data,2,check.column)==2))

  for (i in 1:length(cacos)) {
    caco$cases <- input.data$SampleName[which(input.data[,cacos[i]]==2)]
      cont.ind <- which(input.data[,cacos[i]]==1)
      if (length(caco$cases)*rep>length(cont.ind)) caco$cases <- sample(caco$cases,trunc(length(cont.ind)/rep))
    conts <- sample(cont.ind, length(caco$cases)*rep)
    for (j in 1:rep) {
      caco[[paste0("cont",j)]] <- input.data$SampleName[conts[rep(1:rep,each=length(caco$cases))==j]]
    }
    save(file=paste0(cacos[i],".Rda"), caco)
  }
}

