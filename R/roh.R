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
#'@examples
#'\dontrun{
#'piggyback::pb_download(repo = "soroushmdg/gwid",tag = "v0.0.1",dest = tempdir())
#'ibd_data_file <- paste0(tempdir(),"//chr3.ibd")
#'genome_data_file <- paste0(tempdir(),"//chr3.gds")
#'phase_data_file <- paste0(tempdir(),"//chr3.vcf")
#'case_control_data_file <- paste0(tempdir(),"//case-cont-RA.withmap.Rda")
#'# case-control data
#'case_control <- gwid::case_control(case_control_rda = case_control_data_file)
#'names(case_control) #cases and controls group
#'summary(case_control) # in here, we only consider cases,cont1,cont2,cont3 groups in the study
#'case_control$cases[1:3] # first three subject names of cases group
#'# read SNP data (use SNPRelate to convert it to gds) and count number of minor alleles
#'snp_data_gds <- gwid::build_gwas(gds_data = genome_data_file,
#'caco = case_control,gwas_generator = TRUE)
#'class(snp_data_gds)
#'names(snp_data_gds)
#'head(snp_data_gds$snps) # it has information about counts of minor alleles in each location.
#'# read haplotype data (output of beagle)
#'haplotype_data <- gwid::build_phase(phased_vcf = phase_data_file,caco = case_control)
#'class(haplotype_data)
#'names(haplotype_data)
#'dim(haplotype_data$Hap.1) #22302 SNP and 1911 subjects
#'# read IBD data (output of Refined-IBD)
#'ibd_data <- gwid::build_gwid(ibd_data = ibd_data_file,gwas = snp_data_gds)
#'class(ibd_data)
#'ibd_data$ibd # refined IBD output
#'ibd_data$res # count number of IBD for each SNP location
#'# plot count of IBD in chromosome 3
#'plot(ibd_data,y = c("cases","cont1"),ly = FALSE)
#'# Further investigate location between 117M and 122M
#'# significant number of IBD's in group cases, compare to cont1, cont2 and cont3.
#'plot(ibd_data,y = c("cases","cont1"),snp_start = 117026294,snp_end = 122613594,ly = FALSE)
#'model_fisher <- gwid::fisher_test(ibd_data,case_control,reference = "cases",
#'snp_start = 117026294,snp_end = 122613594)
#'model_permutation <- permutation_test(ibd_data,snp_data_gds,
#'snp_start = 117026294,snp_end = 122613594,nperm=1000,reference = "cases")
#'class(model_fisher)
#'plot(model_fisher, y = c("cases","cont1"),ly = FALSE)
#'hap_str <- gwid::haplotype_structure(ibd_data,phase = haplotype_data,w = 10,
#'snp_start = 117026294,snp_end = 122613594)
#'haplo_freq <- gwid::haplotype_frequency(hap_str)
#'plot(haplo_freq,y = c("cases", "cont1"),plot_type = "haplotype_structure_frequency",
#'nwin = 1, type = "version1",ly = FALSE)
#'}
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
#'@examples
#'\dontrun{
#'piggyback::pb_download(repo = "soroushmdg/gwid",tag = "v0.0.1",dest = tempdir())
#'ibd_data_file <- paste0(tempdir(),"//chr3.ibd")
#'genome_data_file <- paste0(tempdir(),"//chr3.gds")
#'phase_data_file <- paste0(tempdir(),"//chr3.vcf")
#'case_control_data_file <- paste0(tempdir(),"//case-cont-RA.withmap.Rda")
#'# case-control data
#'case_control <- gwid::case_control(case_control_rda = case_control_data_file)
#'names(case_control) #cases and controls group
#'summary(case_control) # in here, we only consider cases,cont1,cont2,cont3 groups in the study
#'case_control$cases[1:3] # first three subject names of cases group
#'# read SNP data (use SNPRelate to convert it to gds) and count number of minor alleles
#'snp_data_gds <- gwid::build_gwas(gds_data = genome_data_file,
#'caco = case_control,gwas_generator = TRUE)
#'class(snp_data_gds)
#'names(snp_data_gds)
#'head(snp_data_gds$snps) # it has information about counts of minor alleles in each location.
#'# read haplotype data (output of beagle)
#'haplotype_data <- gwid::build_phase(phased_vcf = phase_data_file,caco = case_control)
#'class(haplotype_data)
#'names(haplotype_data)
#'dim(haplotype_data$Hap.1) #22302 SNP and 1911 subjects
#'# read IBD data (output of Refined-IBD)
#'ibd_data <- gwid::build_gwid(ibd_data = ibd_data_file,gwas = snp_data_gds)
#'class(ibd_data)
#'ibd_data$ibd # refined IBD output
#'ibd_data$res # count number of IBD for each SNP location
#'# plot count of IBD in chromosome 3
#'plot(ibd_data,y = c("cases","cont1"),ly = FALSE)
#'# Further investigate location between 117M and 122M
#'# significant number of IBD's in group cases, compare to cont1, cont2 and cont3.
#'plot(ibd_data,y = c("cases","cont1"),snp_start = 117026294,snp_end = 122613594,ly = FALSE)
#'model_fisher <- gwid::fisher_test(ibd_data,case_control,reference = "cases",
#'snp_start = 117026294,snp_end = 122613594)
#'model_permutation <- permutation_test(ibd_data,snp_data_gds,
#'snp_start = 117026294,snp_end = 122613594,nperm=1000,reference = "cases")
#'class(model_fisher)
#'plot(model_fisher, y = c("cases","cont1"),ly = FALSE)
#'hap_str <- gwid::haplotype_structure(ibd_data,phase = haplotype_data,w = 10,
#'snp_start = 117026294,snp_end = 122613594)
#'haplo_freq <- gwid::haplotype_frequency(hap_str)
#'plot(haplo_freq,y = c("cases", "cont1"),plot_type = "haplotype_structure_frequency",
#'nwin = 1, type = "version1",ly = FALSE)
#'}
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
#'@examples
#'\dontrun{
#'piggyback::pb_download(repo = "soroushmdg/gwid",tag = "v0.0.1",dest = tempdir())
#'ibd_data_file <- paste0(tempdir(),"//chr3.ibd")
#'genome_data_file <- paste0(tempdir(),"//chr3.gds")
#'phase_data_file <- paste0(tempdir(),"//chr3.vcf")
#'case_control_data_file <- paste0(tempdir(),"//case-cont-RA.withmap.Rda")
#'# case-control data
#'case_control <- gwid::case_control(case_control_rda = case_control_data_file)
#'names(case_control) #cases and controls group
#'summary(case_control) # in here, we only consider cases,cont1,cont2,cont3 groups in the study
#'case_control$cases[1:3] # first three subject names of cases group
#'# read SNP data (use SNPRelate to convert it to gds) and count number of minor alleles
#'snp_data_gds <- gwid::build_gwas(gds_data = genome_data_file,
#'caco = case_control,gwas_generator = TRUE)
#'class(snp_data_gds)
#'names(snp_data_gds)
#'head(snp_data_gds$snps) # it has information about counts of minor alleles in each location.
#'# read haplotype data (output of beagle)
#'haplotype_data <- gwid::build_phase(phased_vcf = phase_data_file,caco = case_control)
#'class(haplotype_data)
#'names(haplotype_data)
#'dim(haplotype_data$Hap.1) #22302 SNP and 1911 subjects
#'# read IBD data (output of Refined-IBD)
#'ibd_data <- gwid::build_gwid(ibd_data = ibd_data_file,gwas = snp_data_gds)
#'class(ibd_data)
#'ibd_data$ibd # refined IBD output
#'ibd_data$res # count number of IBD for each SNP location
#'# plot count of IBD in chromosome 3
#'plot(ibd_data,y = c("cases","cont1"),ly = FALSE)
#'# Further investigate location between 117M and 122M
#'# significant number of IBD's in group cases, compare to cont1, cont2 and cont3.
#'plot(ibd_data,y = c("cases","cont1"),snp_start = 117026294,snp_end = 122613594,ly = FALSE)
#'model_fisher <- gwid::fisher_test(ibd_data,case_control,reference = "cases",
#'snp_start = 117026294,snp_end = 122613594)
#'model_permutation <- permutation_test(ibd_data,snp_data_gds,
#'snp_start = 117026294,snp_end = 122613594,nperm=1000,reference = "cases")
#'class(model_fisher)
#'plot(model_fisher, y = c("cases","cont1"),ly = FALSE)
#'hap_str <- gwid::haplotype_structure(ibd_data,phase = haplotype_data,w = 10,
#'snp_start = 117026294,snp_end = 122613594)
#'haplo_freq <- gwid::haplotype_frequency(hap_str)
#'plot(haplo_freq,y = c("cases", "cont1"),plot_type = "haplotype_structure_frequency",
#'nwin = 1, type = "version1",ly = FALSE)
#'}
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
#'@examples
#'\dontrun{
#'piggyback::pb_download(repo = "soroushmdg/gwid",tag = "v0.0.1",dest = tempdir())
#'ibd_data_file <- paste0(tempdir(),"//chr3.ibd")
#'genome_data_file <- paste0(tempdir(),"//chr3.gds")
#'phase_data_file <- paste0(tempdir(),"//chr3.vcf")
#'case_control_data_file <- paste0(tempdir(),"//case-cont-RA.withmap.Rda")
#'# case-control data
#'case_control <- gwid::case_control(case_control_rda = case_control_data_file)
#'names(case_control) #cases and controls group
#'summary(case_control) # in here, we only consider cases,cont1,cont2,cont3 groups in the study
#'case_control$cases[1:3] # first three subject names of cases group
#'# read SNP data (use SNPRelate to convert it to gds) and count number of minor alleles
#'snp_data_gds <- gwid::build_gwas(gds_data = genome_data_file,
#'caco = case_control,gwas_generator = TRUE)
#'class(snp_data_gds)
#'names(snp_data_gds)
#'head(snp_data_gds$snps) # it has information about counts of minor alleles in each location.
#'# read haplotype data (output of beagle)
#'haplotype_data <- gwid::build_phase(phased_vcf = phase_data_file,caco = case_control)
#'class(haplotype_data)
#'names(haplotype_data)
#'dim(haplotype_data$Hap.1) #22302 SNP and 1911 subjects
#'# read IBD data (output of Refined-IBD)
#'ibd_data <- gwid::build_gwid(ibd_data = ibd_data_file,gwas = snp_data_gds)
#'class(ibd_data)
#'ibd_data$ibd # refined IBD output
#'ibd_data$res # count number of IBD for each SNP location
#'# plot count of IBD in chromosome 3
#'plot(ibd_data,y = c("cases","cont1"),ly = FALSE)
#'# Further investigate location between 117M and 122M
#'# significant number of IBD's in group cases, compare to cont1, cont2 and cont3.
#'plot(ibd_data,y = c("cases","cont1"),snp_start = 117026294,snp_end = 122613594,ly = FALSE)
#'model_fisher <- gwid::fisher_test(ibd_data,case_control,reference = "cases",
#'snp_start = 117026294,snp_end = 122613594)
#'model_permutation <- permutation_test(ibd_data,snp_data_gds,
#'snp_start = 117026294,snp_end = 122613594,nperm=1000,reference = "cases")
#'class(model_fisher)
#'plot(model_fisher, y = c("cases","cont1"),ly = FALSE)
#'hap_str <- gwid::haplotype_structure(ibd_data,phase = haplotype_data,w = 10,
#'snp_start = 117026294,snp_end = 122613594)
#'haplo_freq <- gwid::haplotype_frequency(hap_str)
#'plot(haplo_freq,y = c("cases", "cont1"),plot_type = "haplotype_structure_frequency",
#'nwin = 1, type = "version1",ly = FALSE)
#'}
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
