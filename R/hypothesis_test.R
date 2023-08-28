#' Fisher test
#'
#' @param obj an object
#'
#' @param ... other variables
#'
#' @export
fisher_test <- function(obj, ...) {
  UseMethod("fisher_test")
}

#' Fisher's Exact Test for gwas count data
#'
#' @param obj object of class gwas
#'
#' @param reference reference group of subjects in which we want to perform fisher test
#' test
#' @param snp_start select starting position of snps.
#' @param snp_end select ending position of snp.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#'  You can specify just the initial letter. Only used in the 2 by 2 case
#' @param ... optional arguments to fisher.test
#'
#' @return the output will be a test_snps (data.table) object including 3 columns:
#'  \dQuote{snp_pos}, \dQuote{case_control}, and \dQuote{value} which is
#' a p-values.
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
#'class(model_fisher)
#'plot(model_fisher, y = c("cases","cont1"),ly = FALSE)
#'hap_str <- gwid::haplotype_structure(ibd_data,phase = haplotype_data,w = 10,
#'snp_start = 117026294,snp_end = 122613594)
#'haplo_freq <- gwid::haplotype_frequency(hap_str)
#'plot(haplo_freq,y = c("cases", "cont1"),plot_type = "haplotype_structure_frequency",
#'nwin = 1, type = "version1",ly = FALSE)
#'}
#' @export
fisher_test.gwas <- function(obj, reference, snp_start , snp_end, alternative = c("two.sided","greater","less"), ...) {
  if (missing(obj)){
    stop("please provide gwas object (output of function build_gwas)")
  }

  if (missing(reference)){
    stop("please provide reference names: e.x. 'case1' " )
  }

  alternative <- match.arg(alternative)
  if (sum(unlist(lapply(obj, inherits,"result_snps")))==0){
    obj$snps <- extract(obj)
  }
  snp_pos <- case_control <- NULL
  list_ind_result_snps <- which(unlist(lapply(obj, inherits,"result_snps")))
  snp_pos_total <- unique(obj[[list_ind_result_snps]][["snp_pos"]])
  if (missing(snp_start)) {
    snp_start <- snp_pos_total[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_pos_total[length(snp_pos_total)]
  }
  snp_indx <- which(snp_pos_total >= snp_start & snp_pos_total <= snp_end)
  leni <- diff(range(snp_indx)) + 1
  snp_pos_plot <- snp_pos_total[snp_indx][1:leni]

  snps_mat <- data.table::dcast(obj[[list_ind_result_snps]][snp_pos %in% snp_pos_plot,], snp_pos ~ case_control , value.var = "value")[,-1,with=FALSE]
  snps_mat <- as.matrix(snps_mat)
  nas <- extract(obj,type = "nas")
  nas_mat <- data.table::dcast(nas[snp_pos %in% snp_pos_plot,], snp_pos ~ case_control , value.var = "value")[,-1,with=FALSE]
  nas_mat <- as.matrix(nas_mat)
  idx <- which(colnames(snps_mat) == reference)
  cc <- array(0, dim = c(2, nrow(snps_mat), ncol(snps_mat)))
  pval <- matrix(0, nrow = nrow(snps_mat), ncol = ncol(snps_mat))
  cc[1, , ] <- snps_mat # fill the first matrix with snps_mat
  nums <- as.integer(summary(obj[["caco"]])[1:length(obj[["caco"]])])
  cc[2, , ] <- rep(nums * 2, each = nrow(snps_mat)) - snps_mat - nas_mat * 2
  cc[2,,][cc[2,,] <0] <- 0
  for (i in 1:nrow(snps_mat)) {
    for (j in (1:ncol(snps_mat))[-idx]) pval[i, j] <- stats::fisher.test(cc[, i, c(idx, j)], alternative = alternative, ...)$p.value
  }
  pval[,idx] <- 1
  colnames(pval) <- colnames(snps_mat)
  pval_table <- data.table::as.data.table(pval)[,snp_pos:=snp_pos_total[snp_indx]]
  pval_table <- data.table::melt(pval_table,id.vars="snp_pos",value.name = "value",
                                 measure.vars=colnames(pval),variable.name="case_control")
  output <- pval_table
  class(output) <- append("test_snps",class(output))
  return(output)
}

#' Fisher's Exact Test for gwid count data
#'
#' @param obj An object of class gwid. Output of \code{build_gwid} function
#'
#' @param caco An object of class caco. Output of \code{case_control} function.
#' @param snp_start select starting position of snps.
#' @param snp_end select ending position of snp.
#' @param reference reference group of subjects in which we want to perform fisher test
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#'  You can specify just the initial letter. Only used in the 2 by 2 case
#' @param ... optional arguments to fisher.test
#'
#' @return the output will be a test_snps (data.table) object including 3 columns:
#'  \dQuote{snp_pos}, \dQuote{case_control}, and \dQuote{value} which is
#' a p-values.
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
#'class(model_fisher)
#'plot(model_fisher, y = c("cases","cont1"),ly = FALSE)
#'hap_str <- gwid::haplotype_structure(ibd_data,phase = haplotype_data,w = 10,
#'snp_start = 117026294,snp_end = 122613594)
#'haplo_freq <- gwid::haplotype_frequency(hap_str)
#'plot(haplo_freq,y = c("cases", "cont1"),plot_type = "haplotype_structure_frequency",
#'nwin = 1, type = "version1",ly = FALSE)
#'}
#' @export
fisher_test.gwid <- function(obj, caco,snp_start ,snp_end ,reference , alternative = c("two.sided","greater","less"), ...) {
  if (missing(obj)){
    stop("please provide gwid object (output of function build_gwid)")
  }
  if (missing(caco)){
    stop("please provide caco object(case_control)")
  }
  if (missing(reference)){
    stop("please provide a reference e.x. 'case1' ")
  }
  snp_pos <- NULL
  alternative <- match.arg(alternative)
  #type <- match.arg(type)
  if (sum(unlist(lapply(obj, inherits,"result_snps")))==0){
    obj$res <- obj(obj)
  }
  list_ind_result_snps <- which(unlist(lapply(obj, inherits,"result_snps")))
  snp_pos_total <- unique(obj[[list_ind_result_snps]]$snp_pos)
  if (missing(snp_start)) {
    snp_start <- snp_pos_total[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_pos_total[length(snp_pos_total)]
  }

  snp_indx <- which(snp_pos_total >= snp_start & snp_pos_total <= snp_end)
  leni <- diff(range(snp_indx)) + 1
  snp_pos_plot <- snp_pos_total[snp_indx][1:leni]
  temp <- obj[[list_ind_result_snps]][snp_pos %in% snp_pos_plot,]
  output <- fisher_test(temp,caco=caco,reference=reference,alternative=alternative, ...)
  return(output)
}

#' permutation test
#'
#' @param obj object
#'
#' @param ... other variables
#'
#'
#' @export
permutation_test <- function(obj, ...) {
  UseMethod("permutation_test",obj)
}


#' permutation test for gwid count data
#'
#' @param obj An object of class gwid. Output of \code{build_gwid} function
#'
#' @param gwas object of class gwas
#' @param snp_start select starting position of snps.
#' @param snp_end select ending position of snp.
#' @param nperm Number of permutations.
#' @param reference reference group
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
permutation_test.gwid <- function(obj, gwas, snp_start, snp_end,
                                  nperm = 100, reference = "cases", ...) {
  if (missing(obj)) {
    stop("please provide gwid object (output of function build_gwid)")
  }
  if (missing(gwas)) {
    stop("please provide gwas object")
  }
  if (missing(reference)) {
    stop("please provide reference group")
  }
  if (sum(unlist(lapply(obj, inherits, "result_snps"))) == 0) {
    obj$res <- obj(obj)
  }
  if (missing(snp_start)) {
    snp_start <- snp_pos_total[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_pos_total[length(snp_pos_total)]
  }
  list_ind_result_snps <- which(unlist(lapply(obj, inherits, "result_snps")))
  snp_pos <- case_control <- V1 <- V3 <- NULL
  snp_pos_total <- unique(obj[[list_ind_result_snps]]$snp_pos)
  snp_indx <- which(snp_pos_total >= snp_start & snp_pos_total <= snp_end)
  leni <- diff(range(snp_indx)) + 1
  snp_pos_plot <- snp_pos_total[snp_indx][1:leni]
  temp <- obj[[list_ind_result_snps]][snp_pos %in% snp_pos_plot, ]
  myres <- as.matrix(data.table::dcast(obj$res[snp_pos %in% snp_pos_plot], snp_pos ~ case_control, value.var = "value")[, -1, with = FALSE])
  myres_gap <- myres[,reference] - myres
  perm_matrix <- matrix(0,nrow = nrow(myres_gap),ncol = ncol(myres_gap))
  gwas_temp <- gwas
  gwas_temp$caco <- list(all_subj = unique(unlist(gwas_temp$caco, use.names = FALSE)))
  class(gwas_temp$caco) <- class(gwas[["caco"]])
  build_gwid_modify <- function(gwid = "object of class gwid", gwas = "object of class gwas", gwid_generator = TRUE) {
    ibd <- gwid$ibd
    V1 <- V2 <- V3 <- V4 <- V5 <- V6 <- V7 <- V8 <- V9 <- NULL
    ibd <- ibd[V1 %in% unlist(unique(gwas[["caco"]])) & V3 %in% unlist(unique(gwas[["caco"]]))]
    class(ibd) <- append("IBD", class(ibd))
    seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
    profile <- ind <- vector(mode = "list", length = length(gwas[["caco"]])) # list length 6
    for (j in seq_along(gwas[["caco"]])) {
      ind[[j]] <- which(ibd$V1 %in% gwas[["caco"]][[j]] & ibd$V3 %in% gwas[["caco"]][[j]])
      a1 <- ibd[V1 %in% gwas[["caco"]][[j]] & V3 %in% gwas[["caco"]][[j]]]
      a2 <- seq2(match(a1$V6, gwas$snp.pos), match(a1$V7, gwas$snp.pos))
      Un1 <- unlist(a2)
      profile[[j]] <- Matrix::sparseMatrix(
        i = rep(seq_along(a2), lengths(a2)),
        j = Un1,
        x = 1
      )
      if (ncol(profile[[j]]) < length(gwas$snp.pos)) {
        if (min(Un1) > 1) {
          mytemp <- Matrix(0, nrow = nrow(profile[[j]]), ncol = (min(Un1) - 1))
          profile[[j]] <- cbind(mytemp, profile[[j]])
        }

        if (max(Un1) < length(gwas$snp.pos)) {
          mytemp <- Matrix(0, nrow = nrow(profile[[j]]), ncol = (length(gwas$snp.pos) - max(Un1)))
          profile[[j]] <- cbind(profile[[j]], mytemp)
        }

        profile[[j]] <- methods::as(profile[[j]], "sparseMatrix")
      }
    }

    names(profile) <- names(ind) <- names(gwas[["caco"]])

    class(profile) <- "profile"
    output <- list(profile = profile, IND = ind)
    output$snp_pos <- gwas[["snp.pos"]]
    output$ibd <- ibd
    class(output) <- "gwid"
    if (gwid_generator) {
      output$res <- extract(output)
    }
    return(output)
  }
  myregion2_temp <- build_gwid_modify(gwid = obj, gwas = gwas_temp)
  mres <- myregion2_temp$profile$all_subj
  ibd <- myregion2_temp$ibd
  subj.ind <- list()
  subj.ind[[1]] <- (ibd[,V1])
  subj.ind[[2]] <- (ibd[,V3])
  caco_temp <- gwas[["caco"]]
  mylength <- length(caco_temp[[reference]])
  caco_temp_unlist <- unique(unlist(caco_temp, use.names = FALSE))
  for (i in 1:nperm) {
    cacoi <- matrix(sample(caco_temp_unlist,2*mylength),ncol=2)
    indi <- list()
    mresi <- list()
    indi[[1]] <- which(subj.ind[[1]]%in%cacoi[,1] & subj.ind[[2]]%in%cacoi[,1])
    indi[[2]] <- which(subj.ind[[1]]%in%cacoi[,2] & subj.ind[[2]]%in%cacoi[,2])
    mresi[[1]] <- colSums(mres[indi[[1]],snp_indx])
    mresi[[2]] <- colSums(mres[indi[[2]],snp_indx])
    perm_matrix <- perm_matrix + 1*( myres_gap - matrix(mresi[[1]]-mresi[[2]],length(mresi[[1]]-mresi[[2]]),ncol(myres_gap))<= 0 )
  }
  perm_matrix <- cbind(snp_pos_plot,perm_matrix/nperm)
  colnames(perm_matrix) <- c("snp_pos",colnames(perm_matrix)[-1])
  perm_table <- data.table::as.data.table(perm_matrix)
  pval <- data.table::melt(perm_table, id.vars = "snp_pos", variable.name = "case_control", value.name = "value")
  class(pval) <- append("test_snps", class(pval))
  return(pval)
}

#' Permutation test for gwas object
#'
#' @param obj object of class gwas
#'
#' @param snp_start elect starting position of snps.
#' @param snp_end select ending position of snp.
#' @param nperm Number of permutations.
#' @param reference reference group of subjects in which we want to perform fisher test
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
permutation_test.gwas <- function(obj, snp_start, snp_end,
                                  nperm = 1000, reference = "cases", ...) {
  if (missing(obj)) {
    stop("please provide gwas object")
  }

  if (missing(snp_start)) {
    snp_start <- snp_pos_total[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_pos_total[length(snp_pos_total)]
  }
  snp_pos <- case_control <- NULL
  list_ind_result_snps <- which(unlist(lapply(obj, inherits, "result_snps")))
  snp_pos_total <- unique(obj[[list_ind_result_snps]]$snp_pos)
  snp_indx <- which(snp_pos_total >= snp_start & snp_pos_total <= snp_end)
  leni <- diff(range(snp_indx)) + 1
  snp_pos_plot <- snp_pos_total[snp_indx][1:leni]
  caco <- obj$caco
  total_subjects <- sapply(obj$caco, length)
  myres <- as.matrix(data.table::dcast(obj$snps[snp_pos %in% snp_pos_plot], snp_pos ~ case_control, value.var = "value")[, -1, with = FALSE])
  myres_gap <- myres[, reference] - myres
  perm_matrix <- matrix(0, nrow = nrow(myres_gap), ncol = ncol(myres_gap))
  gwas_temp <- obj
  gwas_temp$caco <- list(all_subj = (unlist(gwas_temp$caco, use.names = FALSE)))
  class(gwas_temp$caco) <- class(obj$caco)
  mres <- do.call(what = "rbind", args = obj$smp.snp)
  gwas_temp$smp.snp <- list(all_subj = mres)
  for (i in 1:nperm) {
    caco_temp <- obj$caco
    cacoi <- split(sample(gwas_temp$caco$all_subj), rep(1:length(obj$caco), total_subjects))
    names(cacoi) <- names(obj$caco)
    new_snps <- sapply(cacoi, function(x) {
      colSums(gwas_temp$smp.snp$all_subj[rownames(gwas_temp$smp.snp$all_subj) %in% x, snp_indx])
    })
    myres_gap_perm <- new_snps[, reference] - new_snps
    perm_matrix <- perm_matrix + 1 * (myres_gap - myres_gap_perm <= 0)
  }
  perm_matrix <- cbind(snp_pos_plot, perm_matrix / nperm)
  colnames(perm_matrix) <- c("snp_pos", colnames(perm_matrix)[-1])
  perm_table <- data.table::as.data.table(perm_matrix)
  pval <- data.table::melt(perm_table, id.vars = "snp_pos", variable.name = "case_control", value.name = "value")
  class(pval) <- append("test_snps", class(pval))
  return(pval)
}

G_Test <- function (x, y = NULL, correct = c("none", "williams", "yates"),
                    p = rep(1/length(x), length(x)))
{
  DNAME <- deparse(substitute(x))
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.matrix(x)) {
    if (min(dim(x)) == 1)
      x <- as.vector(x)
  }
  if (!is.matrix(x) && !is.null(y)) {
    if (length(x) != length(y))
      stop("x and y must have the same length")
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    OK <- stats::complete.cases(x, y)
    x <- as.factor(x[OK])
    y <- as.factor(y[OK])
    if ((nlevels(x) < 2) || (nlevels(y) < 2))
      stop("x and y must have at least 2 levels")
    x <- table(x, y)
  }
  if (any(x < 0) || any(is.na(x)))
    stop("all entries of x must be nonnegative and finite")
  if ((n <- sum(x)) == 0)
    stop("at least one entry of x must be positive")
  correct <- match.arg(correct)
  if (is.matrix(x)) {
    nrows <- nrow(x)
    ncols <- ncol(x)
    if (correct == "yates") {
      if (dim(x)[1] != 2 || dim(x)[2] != 2)
        stop("Yates' correction requires a 2 x 2 matrix")
      if ((x[1, 1] * x[2, 2]) - (x[1, 2] * x[2, 1]) > 0) {
        x <- x + 0.5
        diag(x) <- diag(x) - 1
      }
      else {
        x <- x - 0.5
        diag(x) <- diag(x) + 1
      }
    }
    sr <- apply(x, 1, sum)
    sc <- apply(x, 2, sum)
    E <- outer(sr, sc, "*")/n
    g <- 0
    for (i in 1:nrows) {
      for (j in 1:ncols) {
        if (x[i, j] != 0)
          g <- g + x[i, j] * log(x[i, j]/E[i, j])
      }
      q <- 1
      if (correct == "williams") {
        row.tot <- col.tot <- 0
        for (i in 1:nrows) {
          row.tot <- row.tot + 1/(sum(x[i, ]))
        }
        for (j in 1:ncols) {
          col.tot <- col.tot + 1/(sum(x[, j]))
        }
        q <- 1 + ((n * row.tot - 1) * (n * col.tot -
                                         1))/(6 * n * (ncols - 1) * (nrows - 1))
      }
      STATISTIC <- G <- 2 * g/q
      PARAMETER <- (nrow(x) - 1) * (ncol(x) - 1)
      PVAL <- 1 - stats::pchisq(STATISTIC, df = PARAMETER)
      if (correct == "none")
        METHOD <- "Log likelihood ratio (G-test) test of independence without correction"
      if (correct == "williams")
        METHOD <- "Log likelihood ratio (G-test) test of independence with Williams' correction"
      if (correct == "yates")
        METHOD <- "Log likelihood ratio (G-test) test of independence with Yates' correction"
    }
  }
  else {
    METHOD <- "Log likelihood ratio (G-test) goodness of fit test"
    if (length(x) == 1)
      stop("x must at least have 2 elements")
    if (length(x) != length(p))
      stop("x and p must have the same number of elements")
    E <- n * p
    if (correct == "yates") {
      if (length(x) != 2)
        stop("Yates' correction requires 2 data values")
      if ((x[1] - E[1]) > 0.25) {
        x[1] <- x[1] - 0.5
        x[2] <- x[2] + 0.5
      }
      else if ((E[1] - x[1]) > 0.25) {
        x[1] <- x[1] + 0.5
        x[2] <- x[2] - 0.5
      }
    }
    names(E) <- names(x)
    g <- 0
    for (i in 1:length(x)) {
      if (x[i] != 0)
        g <- g + x[i] * log(x[i]/E[i])
    }
    q <- 1
    if (correct == "williams") {
      q <- 1 + (length(x) + 1)/(6 * n)
    }
    STATISTIC <- G <- 2 * g/q
    PARAMETER <- length(x) - 1
    PVAL <- stats::pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  }
  names(STATISTIC) <- "G"
  names(PARAMETER) <- "X-squared df"
  names(PVAL) <- "p.value"
  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, method = METHOD, data.name = DNAME, observed = x,
                 expected = E), class = "htest")
}

#' perform gtest
#'
#' @param haplotype_structure object of a class
#'
#' @param ... other variables
#'
#' @export
gtest <- function(haplotype_structure, ...) {
  UseMethod("gtest")
}

#' Perform G-test on haplotype structures extracted from \code{haplotype_structure} function
#'
#' @param haplotype_structure An object of class haplotype_structure. Output of
#' \code{haplotype_structure} function.
#'
#' @param reference reference group of subjects in which we want to perform G-test
#'
#' @param ... other variables
#'
#'@return the output will be a test_snps (data.table) object including 3 columns:
#'  \dQuote{snp_pos}, \dQuote{case_control}, and \dQuote{value} which is
#' a p-values.
#'
#' @export
gtest.haplotype_structure <- function(haplotype_structure, reference, ...) {

  if (missing(haplotype_structure) | missing(reference)){
    stop("please provide function arguments")
  }
  pval <- matrix(0, nrow = length(unique(haplotype_structure$snp_pos)), ncol = length(unique(haplotype_structure$case_control)))
  colnames(pval) <- levels(haplotype_structure$case_control)
  idx <- which(colnames(pval) == reference)
  for (j in (1:ncol(pval))[-idx]) {
    caco <- c(reference, colnames(pval)[j])
    case_control <- snp_pos <- structures <- window_number <- NULL
    tmp <- haplotype_structure[case_control %in% caco, list(case_control, structures, window_number)]

    tmp <- split(tmp, tmp[, window_number])

    tmp <- lapply(tmp, function(x) {
      y <- x[case_control %in% caco, !"window_number"]
      y <- table(y)[caco, ]
    })
    pval[, j] <- unlist(lapply(tmp, function(x) {
      y <- G_Test(x, correct = "williams")$p.value
      y
    }))
  }
  pval[, idx] <- 1
  pval_table <- data.table::as.data.table(pval)[, snp_pos := sort(unique(haplotype_structure$snp_pos))]
  pval_table <- data.table::melt(pval_table,
                                 id.vars = "snp_pos", value.name = "value",
                                 measure.vars = colnames(pval), variable.name = "case_control"
  )
  output <- data.table::copy(pval_table)
  class(output) <- append("test_snps",class(output))
  return(output)
}


#' fisher exact test for result_snps count data
#'
#' @param obj An object of class result_snps
#'
#' @param caco An object of class caco. Output of \code{case_control} function.
#' @param reference reference group of subjects in which we want to perform fisher test.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#'  You can specify just the initial letter. Only used in the 2 by 2 case
#' @param ... optional arguments to fisher.test
#' @return the output will be a test_snps (data.table) object including 3 columns:
#'  \dQuote{snp_pos}, \dQuote{case_control}, and \dQuote{value} which is
#' a p-values.
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
fisher_test.result_snps <- function(obj, caco, reference, alternative = c("two.sided","greater","less"), ...) {

  if (missing(obj)){
    stop("Please provide an object of class result_snps")
  }
  if (missing(caco)){
    stop("Please provide an object of class caco")
  }
  if (missing(reference)){
    stop("Please provide reference e.x. 'case1' ")
  }
  alternative <- match.arg(alternative)
  res_mat <- data.table::dcast(obj, snp_pos ~ case_control , value.var = "value")
  snp_pos <- res_mat[["snp_pos"]]
  res_mat <- as.matrix(res_mat[,-1,with=FALSE])

  #idx is column index of reference
  idx <- which(colnames(res_mat) == reference)
  cc <- array(0, dim = c(2, nrow(res_mat), ncol(res_mat)))
  pval <- matrix(0, nrow = nrow(res_mat), ncol = ncol(res_mat))
  cc[1, , ] <- res_mat
  nums <- as.integer(summary(caco)[1:length(caco)]) # number of samples in each case_control

  cc[2, , ] <- rep(nums * (nums - 1) / 2, each = nrow(res_mat)) - res_mat

  for (i in 1:nrow(res_mat)) {
    for (j in (1:ncol(res_mat))[-idx]) pval[i, j] <- fisher.test(cc[, i, c(idx, j)], alternative = alternative)$p.value
  }
  colnames(pval) <- colnames(res_mat)
  pval[,idx] <- 1
  pval_table <- data.table::as.data.table(pval)[,snp_pos:=snp_pos]
  # transform pval to long format 3 columns (snp_pos, case_control, value)
  pval_table <- data.table::melt(pval_table,id.vars="snp_pos",value.name = "value",
                                 measure.vars=colnames(pval),variable.name="case_control")
  class(pval_table) <- append("test_snps",class(pval_table))
  return(pval_table)
}



#' Permutation test for `haplotype_structure` object
#'
#' @param obj object of class `haplotype_structure`
#'
#' @param nperm Number of permutations.
#' @param reference reference group of subjects in which we want to perform `gtest`
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
permutation_test.haplotype_structure <- function(obj, nperm, reference, ...) {
  if (missing(obj) | missing(reference)) {
    stop("please provide function arguments")
  }
  statistics <- matrix(0, nrow = length(unique(obj$snp_pos)), ncol = length(unique(obj$case_control)))
  colnames(statistics) <- levels(obj$case_control)
  idx <- which(colnames(statistics) == reference)
  statistics_perm <- vector(mode = "list", length = nperm)
  for (i in 1:nperm) {
    statistics_perm[[i]] <- matrix(0, nrow = nrow(statistics), ncol = ncol(statistics))
  }
  Freq <- NULL
  for (j in (1:ncol(statistics))[-idx]) {
    caco <- c(reference, colnames(statistics)[j])
    case_control <- snp_pos <- structures <- window_number <- NULL
    tmp <- obj[case_control %in% caco, list(case_control, structures, window_number)]

    tmp <- split(tmp, tmp[, window_number])

    tmp <- lapply(tmp, function(x) {
      y <- x[case_control %in% caco, !"window_number"]
      y <- table(y)[caco, ]
    })
    statistics[, j] <- unlist(lapply(tmp, function(x) {
      y <- G_Test(x, correct = "williams")$statistic
      y
    }))

    for (i in 1:nperm) {
      statistics_perm[[i]][, j] <- unlist(lapply(tmp, function(x) {
        data_frame_tmp <- as.data.frame(x)
        data_frame_tmp[, 3] <- sample(data_frame_tmp[, 3])
        cn <- colnames(data_frame_tmp)
        shuffled_table_obj <- xtabs(Freq ~ ., data_frame_tmp)
        y <- G_Test(shuffled_table_obj, correct = "williams")$statistic
        y
      }))
    }
  }

  pval <- matrix(0, nrow = nrow(statistics), ncol = ncol(statistics))
  pval <- Reduce("+", lapply(statistics_perm, function(x) {
    x > statistics
  })) / nperm
  colnames(pval) <- colnames(statistics)
  pval[, idx] <- 1
  pval_table <- data.table::as.data.table(pval)[, snp_pos := sort(unique(obj$snp_pos))]
  pval_table <- data.table::melt(pval_table,
    id.vars = "snp_pos", value.name = "value",
    measure.vars = colnames(statistics), variable.name = "case_control"
  )
  output <- data.table::copy(pval_table)
  class(output) <- append("test_snps", class(output))
  return(output)
}
