#' haplotype structures in a window
#'
#' @param obj object
#'
#' @param ... other variables
#' @return The output will be an object of class haplotype_structure (data.table) that has information
#' about subjects haplotype structures in a a window.
#' @export
haplotype_structure <- function(obj, ...) {
  UseMethod("haplotype_structure")
}

#' extract haplotype structures of pairwise ibd samples in a window
#'
#' @param obj An object of class gwid. Output of \code{build_gwid} function.
#'
#' @param phase An object of class phase. Output of \code{build_phase} function.
#' @param w window size
#' @param snp_start select starting position of snps.
#' @param snp_end select ending position of snps.
#' @param ... other variables
#'
#' @return The output will be an object of class haplotype_structure (data.table) that has information
#' about subjects haplotype structures in a a window.
#'
#'@examples
#'\donttest{
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
#'plot(ibd_data,y = c("cases","cont1"),snp_start = 119026294,snp_end = 120613594,ly = FALSE)
#'model_fisher <- gwid::fisher_test(ibd_data,case_control,reference = "cases",
#'snp_start = 119026294,snp_end = 120613594)
#'class(model_fisher)
#'plot(model_fisher, y = c("cases","cont1"),ly = FALSE)
#'hap_str <- gwid::haplotype_structure(ibd_data,phase = haplotype_data,w = 10,
#'snp_start = 119026294,snp_end = 120613594)
#'haplo_freq <- gwid::haplotype_frequency(hap_str)
#'plot(haplo_freq,y = c("cases", "cont1"),plot_type = "haplotype_structure_frequency",
#'nwin = 1, type = "version1",ly = FALSE)
#'}
#' @export
haplotype_structure.gwid <- function(obj, phase, w = 10, snp_start, snp_end, ...) {
  #browser()
  if (missing(obj)) {
    stop("please provide gwid object (output of function build_gwid)")
  }
  snp_pos_total <- unique(obj[[
    which(unlist(lapply(obj, inherits, "result_snps")))
  ]]$snp_pos)
  if (missing(snp_start)) {
    snp_start <- snp_pos_total[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_pos_total[length(snp_pos_total)]
  }
  smp_index <- window_number <- snp_indx <- haplo_index <-
    snp_pos <- structures <- smp <- snp_index <- sum_window <- NULL

  snp_indx <- which(snp_pos_total >= snp_start & snp_pos_total <= snp_end)
  leni <- diff(range(snp_indx)) - w + 2
  snp_pos_plot <- snp_pos_total[snp_indx][1:leni]
  snp_pos_plot <- data.table::data.table(snp_pos = snp_pos_plot)
  snp_pos_plot <- cbind(snp_pos_plot, window_number = seq(1, nrow(snp_pos_plot)))
  if (sum(snp_indx)==0){
    stop("length of select region should be positive")
  }

  if (w >= diff(range(snp_indx))) {
    stop("window size should be smaller than number of snps")
  }

  if (length(snp_indx) == length(snp_pos_total)) {
    Mres_reduced <- lapply(obj[[
      which(unlist(lapply(obj, inherits, "profile")))
    ]], as.matrix)
  } else {
    Mres_reduced <- lapply(obj[[
      which(unlist(lapply(obj, inherits, "profile")))
    ]], "[", j = snp_indx,exact=TRUE)
  }
  df <- lapply(Mres_reduced, function(x) {
    y <- t(apply(x, 1, RcppRoll::roll_sum, n = w))
    y <- data.table::as.data.table(y)
    return(y)
  })

  for (i in 1:length(df)) {
    haplotype <- obj[[
      which(unlist(lapply(obj, inherits, "IBD")))
    ]][obj$IND[[i]], c(2, 1)]
    data.table::setnames(haplotype, c("V2", "V1"), c("haplo_index", "smp"))
    df[[i]] <- cbind(haplotype, df[[i]])
  }

  sample_indexing <- data.table::as.data.table(x = colnames(phase[[1]]))[, smp_index := .I]
  data.table::setnames(sample_indexing, old = "V1", new = "smp")
  data.table::setcolorder(sample_indexing, c("smp_index", setdiff(names(sample_indexing), "smp_index")))
  data.table::setkey(sample_indexing, smp)

  df <- lapply(df, function(x) {
    y <- data.table::melt(x, id.vars = c("haplo_index", "smp"), variable.name = "window_number", value.name = "sum_window")
    y <- y[sum_window == w][, window_number := lapply(.SD, function(x) y <- as.integer(substring(x, 2))), .SDcols = "window_number"][, sum_window := NULL, ]
    return(y)
  })

  df <- data.table::rbindlist(l = df, use.names = TRUE, idcol = "case_control")
  data.table::setkey(df, smp)

  df <- data.table::merge.data.table(x = df, y = sample_indexing, by = "smp")

  data.table::setkey(df, smp_index, window_number)
  df[, snp_index := window_number + min(snp_indx) - min(window_number)]
  data.table::setcolorder(df, c("case_control", "haplo_index", "smp", "smp_index", "snp_index", "window_number"))

  result <- Matrix::Matrix(0, nr = nrow(df), nc = w)
  result1 <- data.table::as.data.table(matrix(0, nrow = nrow(df), ncol = w))
  my_phase1 <- phase[[1]][range(snp_indx)[1]:nrow(phase[[1]]), ]
  my_phase2 <- phase[[2]][range(snp_indx)[1]:nrow(phase[[2]]), ]


  df1 <- data.table::copy(df)
  for (i in 1:w) {
    tmp1 <- data.table::as.data.table(cbind(which(df$haplo_index == 1), my_phase1[as.matrix(df1[haplo_index == 1, c(6, 4)])]))
    tmp2 <- data.table::as.data.table(cbind(which(df$haplo_index == 2), my_phase2[as.matrix(df1[haplo_index == 2, c(6, 4)])]))
    data.table::setnames(x = tmp1, old = c("V1", "V2"), new = c("ind", "val"))
    data.table::setnames(x = tmp2, old = c("V1", "V2"), new = c("ind", "val"))
    data.table::setkey(tmp1)
    data.table::setkey(tmp2)

    result1[, i] <- data.table::merge.data.table(x = tmp1, y = tmp2, by = c("ind", "val"), all = TRUE)[, 2]
    df1[, 6] <- df1[, 6] + 1
  }
  rm(df1)

  res <- cbind(df, structures = result1[, apply(.SD, 1, paste, collapse = ""), .SDcols = colnames(result1)])
  res <- res[, case_control := factor(case_control, levels = names(obj[[
    which(unlist(lapply(obj, inherits, "profile")))
  ]]))]
  data.table::setkey(res, window_number)
  res <- merge(res, snp_pos_plot, by = "window_number")
  res <- res[, list(case_control, snp_pos, window_number, smp, structures)]
  output <- res
  class(output) <- append("haplotype_structure", class(output))
  return(output)
}


#' haplotype frequency
#'
#' @param haplotype_structure object of class haplotype structure
#'
#' @param ... other variables
#' @return An object of class haplotype_frequency contains of two objects. first one
#' is object of haplotype_structure_frequency (data.table) and second one is object of class result_snps(data.table)
#' @examples
#'\donttest{
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
#'plot(ibd_data,y = c("cases","cont1"),snp_start = 119026294,snp_end = 120613594,ly = FALSE)
#'model_fisher <- gwid::fisher_test(ibd_data,case_control,reference = "cases",
#'snp_start = 119026294,snp_end = 120613594)
#'class(model_fisher)
#'plot(model_fisher, y = c("cases","cont1"),ly = FALSE)
#'hap_str <- gwid::haplotype_structure(ibd_data,phase = haplotype_data,w = 10,
#'snp_start = 119026294,snp_end = 120613594)
#'haplo_freq <- gwid::haplotype_frequency(hap_str)
#'plot(haplo_freq,y = c("cases", "cont1"),plot_type = "haplotype_structure_frequency",
#'nwin = 1, type = "version1",ly = FALSE)
#'}
#' @export
haplotype_frequency <- function(haplotype_structure, ...) {
  UseMethod("haplotype_frequency")
}

#' haplotype frequency in sliding windows
#'
#' @param haplotype_structure An object of class haplotype_structure. Output of
#' \code{haplotype_structure} function.
#'
#' @param ... other variables
#'
#' @return An object of class haplotype_frequency contains of two objects. first one
#' is object of haplotype_structure_frequency (data.table) and second one is object of class result_snps(data.table)
#'
#'@examples
#'\donttest{
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
#'plot(ibd_data,y = c("cases","cont1"),snp_start = 119026294,snp_end = 120613594,ly = FALSE)
#'model_fisher <- gwid::fisher_test(ibd_data,case_control,reference = "cases",
#'snp_start = 119026294,snp_end = 120613594)
#'class(model_fisher)
#'plot(model_fisher, y = c("cases","cont1"),ly = FALSE)
#'hap_str <- gwid::haplotype_structure(ibd_data,phase = haplotype_data,w = 10,
#'snp_start = 119026294,snp_end = 120613594)
#'haplo_freq <- gwid::haplotype_frequency(hap_str)
#'plot(haplo_freq,y = c("cases", "cont1"),plot_type = "haplotype_structure_frequency",
#'nwin = 1, type = "version1",ly = FALSE)
#'}
#' @export
haplotype_frequency.haplotype_structure <- function(haplotype_structure, ...) {
  if (missing(haplotype_structure)) {
    stop("please provide object of class haplotype_structure")
  }
  window_number <- structures <- case_control <- NULL
  df <- haplotype_structure[, .N, by = list(window_number, case_control, structures)][, window_number := factor(window_number, levels = unique(window_number))]
  data.table::setnames(df, "N", "n")
  df <- split(df, f = df[, window_number])
  df <- lapply(df, function(x) {
    x[, structures := factor(structures, levels = unique(structures))]
    x
  })
  tmp <- data.table::rbindlist(l = df, use.names = TRUE)

  df <- lapply(df, function(x) {
    y <- x[
      data.table::CJ(window_number, case_control, structures, unique = T),
      on = list(window_number, case_control, structures)
    ]

    data.table::setnafill(y, fill = 0, cols = "n")
  })
  df <- data.table::rbindlist(l = df, use.names = TRUE)
  df[, window_number := as.integer(window_number)]
  class(df) <- append("haplotype_structure_frequency", class(df))
  n <- snp_pos <- value <- NULL
  df1 <- tmp[, list(value = max(n)), by = list(window_number, case_control)]
  df1 <- split(df1, f = df1[, case_control])
  df1 <- lapply(df1, function(x) {
    cbind(x, snp_pos = sort(unique(haplotype_structure$snp_pos)))
  })
  df1 <- data.table::rbindlist(l = df1, use.names = TRUE)
  df1 <- df1[, list(snp_pos, case_control, value)]
  class(df1) <- append("result_snps", class(df1))
  df2 <- tmp[, list(value = sum(n)), by = list(window_number, case_control)]
  df2 <- split(df2, f = df2[, case_control])
  df2 <- lapply(df2, function(x) {
    cbind(x, snp_pos = sort(unique(haplotype_structure$snp_pos)))
  })
  df2 <- data.table::rbindlist(l = df2, use.names = TRUE)
  df2 <- df2[, list(snp_pos, case_control, value)]
  class(df2) <- append("result_snps", class(df2))


  output <- list(hap_freq = df, hap_most_freq = df1)
  class(output) <- append("haplotype_frequency", class(output))
  return(output)
}

#' extract haplotype structures of individuals in a window
#'
#' @param obj object of class gwas
#'
#' @param phase An object of class phase. Output of \code{build_phase} function
#' @param w window size
#' @param snp_start select starting position of snps.
#' @param snp_end select ending position of snps.
#' @param ... other variables
#'
#' @return The output will be an object of class haplotype_structure (data.table) that has information
#' about subjects haplotype structures in a a window.
#'
#'@examples
#'\donttest{
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
#'plot(ibd_data,y = c("cases","cont1"),snp_start = 119026294,snp_end = 120613594,ly = FALSE)
#'model_fisher <- gwid::fisher_test(ibd_data,case_control,reference = "cases",
#'snp_start = 119026294,snp_end = 120613594)
#'class(model_fisher)
#'plot(model_fisher, y = c("cases","cont1"),ly = FALSE)
#'hap_str <- gwid::haplotype_structure(ibd_data,phase = haplotype_data,w = 10,
#'snp_start = 119026294,snp_end = 120613594)
#'haplo_freq <- gwid::haplotype_frequency(hap_str)
#'plot(haplo_freq,y = c("cases", "cont1"),plot_type = "haplotype_structure_frequency",
#'nwin = 1, type = "version1",ly = FALSE)
#'}
#' @export
haplotype_structure.gwas <- function(obj, phase, w = 10, snp_start, snp_end, ...) {
  if (missing(phase)) {
    stop("please provide phase object")
  }
  if (missing(snp_start)) {
    snp_start <- obj$snp.pos[1]
  }
  if (missing(snp_end)) {
    snp_end <- obj$snp.pos[length(obj$snp.pos)]
  }
  snp_index <- which(obj$snp.pos >= snp_start & obj$snp.pos <= snp_end)
  leni <- diff(range(snp_index)) - w + 2
  snp_pos_plot <- obj$snp.pos[snp_index][1:leni]
  if (w >= diff(range(snp_index))) {
    stop("window size should be smaller than number of snps")
  }
  if (w >= length(snp_index)) {
    stop("window size should be smaller than number of snps")
  }
  smp_index <- window_number <- snp_pos <- case_control <- smp <- .id <- structures <- n <- NULL
  sample_indexing <- data.table::as.data.table(x = colnames(phase[[1]]))[, smp_index := .I]
  data.table::setnames(sample_indexing, old = "V1", new = "smp")
  data.table::setcolorder(sample_indexing, c("smp_index", setdiff(names(sample_indexing), "smp_index")))
  my_phase1 <- phase[[1]][range(snp_index)[1]:nrow(phase[[1]]), ]
  my_phase2 <- phase[[2]][range(snp_index)[1]:nrow(phase[[2]]), ]
  sample_index_caco <- list()
  for (i in 1:length(obj[[
    which(unlist(lapply(obj, inherits, "caco")))
  ]])) {
    snp_index_dt <- data.table::as.data.table(snp_index)[1:leni][, window_number := snp_index - min(snp_index) + 1][, snp_pos := snp_pos_plot]
    sample_index_caco[[i]] <- sample_indexing[smp %in% obj[[
      which(unlist(lapply(obj, inherits, "caco")))
    ]][[i]], ]
    n <- nrow(sample_index_caco[[i]])
    sample_index_caco[[i]] <- cbind(sample_index_caco[[i]], snp_index_dt[rep(seq_len(length.out = dim(snp_index_dt)[1]), n), ])
  }
  names(sample_index_caco) <- names(obj[[
    which(unlist(lapply(obj, inherits, "caco")))
  ]])

  .id <- NULL
  myind <- data.table::rbindlist(sample_index_caco, idcol = TRUE)[, case_control := .id][, .id := NULL]

  result1 <- Matrix::Matrix(NA, nr = nrow(myind), nc = w)
  result2 <- Matrix::Matrix(NA, nr = nrow(myind), nc = w)

  for (i in 1:w) {
    result1[, i] <- my_phase1[as.matrix(myind[, c(4, 1)])]
    result2[, i] <- my_phase2[as.matrix(myind[, c(4, 1)])]
    myind[, 4] <- myind[, 4] + 1
  }

  myind[, 4] <- myind[, 4] - w
  res1 <- cbind(myind, apply(result1 * 1, 1, paste, collapse = ""))
  res2 <- cbind(myind, apply(result2 * 1, 1, paste, collapse = ""))
  res <- rbind(res1, res2)
  data.table::setnames(res, old = "V2", new = "structures")
  res <- res[, case_control := factor(case_control, levels = names(obj$caco))]
  data.table::setcolorder(res, c("case_control", "snp_pos", "window_number", "smp", "structures", "smp_index", "snp_index"))
  # browser()
  res <- res[, list(case_control, snp_pos, window_number, smp, structures)]

  class(res) <- append("haplotype_structure", class(res))
  return(res)
}
