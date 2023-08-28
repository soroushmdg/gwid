#' Extract information from SNP GDS file.
#'
#' @param obj an object of class gwas
#' @param ... other arguments
#'
#' @return extract object instants
#' @export
extract <- function(obj, ...) {
  UseMethod("extract")
}

#' Extract information from SNP GDS file.
#'
#' @param obj object of class gwas.
#'
#' @param type indicate type of aggregation on sample-snp data and must be
#' one of snps, snp2, or nas
#'
#' @param snp_start select starting position of snp, which we want to aggregate.
#' @param snp_end select ending position of snp, which we want to aggregate.
#' @param ... other arguments
#'
#' @return the output will be a result_snps (data.table) object including 3 columns
#' including, snp_pos, case_control, and value
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
extract.gwas <- function(obj, type = c("snps", "snp2", "nas"), snp_start, snp_end, ...) {
  # summation of each column of smp_snp matrix result in 5296*6 matrix
  type <- match.arg(type)
  output <- list()
  if ("snps" %in% type) {
    output[["snps"]] <- matrix(unlist(lapply(obj[["smp.snp"]], colSums, na.rm = TRUE)),
      ncol = length(obj[["caco"]])
    )
  }
  if ("snp2" %in% type) {
    output[["snp2"]] <- matrix(unlist(lapply(lapply(obj[["smp.snp"]], function(x) x == 2), colSums,
      na.rm = TRUE
    )), ncol = length(obj[["caco"]]))
  }
  if ("nas" %in% type) {
    output[["nas"]] <- matrix(unlist(lapply(lapply(obj[["smp.snp"]], function(x) is.na(x)), colSums)),
      ncol = length(obj[["caco"]])
    )
  }
  output <- lapply(output, function(x) {
    mode(x) <- "integer" # kind of like type (storage mode)
    colnames(x) <- names(obj[["caco"]])
    x
  })
  output <- lapply(output, data.table::as.data.table)
  output <- Map(cbind, output, snp_pos = list(obj[["snp.pos"]]))
  output <- lapply(output, data.table::melt, id.vars = "snp_pos", variable.name = "case_control", value.name = "value")
  snp_pos <- case_control <- value <- NULL
  output <- lapply(output, function(x) {
    snp_pos <- colnames(x)[2]
    data.table::setkey(x, snp_pos)
    class(x) <- append("result_snps", class(x))
    return(x)
  })
  output2 <- output[[type]]
  if (!missing(snp_start)) {
    output2 <- output2[snp_pos >= snp_start]
  }
  if (!missing(snp_end)) {
    output2 <- output2[snp_pos <= snp_end]
  }
  return(output2)
}

#' Extract information from ibd data.
#'
#' @param obj object of class gwid(output of function build_gwid)
#'
#' @param snp_start select starting position of snp, which we want to aggregate.
#' @param snp_end select ending position of snp, which we want to aggregate.
#' @param ... other objects
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
#'class(model_fisher)
#'plot(model_fisher, y = c("cases","cont1"),ly = FALSE)
#'hap_str <- gwid::haplotype_structure(ibd_data,phase = haplotype_data,w = 10,
#'snp_start = 117026294,snp_end = 122613594)
#'haplo_freq <- gwid::haplotype_frequency(hap_str)
#'plot(haplo_freq,y = c("cases", "cont1"),plot_type = "haplotype_structure_frequency",
#'nwin = 1, type = "version1",ly = FALSE)
#'}
#' @export
extract.gwid <- function(obj = "object of class gwid", snp_start, snp_end, ...) {
  list_ind_profile <- which(unlist(lapply(obj, inherits, "profile")))
  res <- sapply(obj[[list_ind_profile]], colSums)
  res <- cbind(snp_pos = obj[["snp_pos"]], res)
  res <- data.table::as.data.table(x = res)
  snp_pos <- NULL
  data.table::setkey(x = res, snp_pos)
  res <- data.table::melt(res, id.vars = "snp_pos", variable.name = "case_control", value.name = "value")
  if (!missing(snp_start)) {
    res <- res[snp_pos >= snp_start]
  }
  if (!missing(snp_end)) {
    res <- res[snp_pos <= snp_end]
  }
  class(res) <- append("result_snps", class(res))
  return(res)
}


#' extract component of an object
#'
#' @param obj obj
#'
#' @param ... other variables
#'
#' @export
extract_window <- function(obj, ...) {
  UseMethod("extract_window")
}


#' Extract information from ibd data in a moving window
#'
#' @param obj object of class gwid(output of function build_gwid)
#'
#' @param w window size
#' @param snp_start select starting position of snp, which we want to aggregate.
#' @param snp_end select ending position of snp, which we want to aggregate.
#' @param ... other variables
#' @return the output will be a result_snps (data.table) object including 3 columns
#' including, \dQuote{snp_pos}, \dQuote{case_control}, and \dQuote{value}
#'
#' @examples
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
#'
#' @export
extract_window.gwid <- function(obj, w = 10, snp_start, snp_end, ...) {
  snp_pos_total <- unique(obj[[
    which(unlist(lapply(obj, inherits, "result_snps")))
  ]][["snp_pos"]])
  if (missing(snp_start)) {
    snp_start <- snp_pos_total[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_pos_total[length(snp_pos_total)]
  }

  snp_indx <- which(snp_pos_total >= snp_start & snp_pos_total <= snp_end)
  leni <- diff(range(snp_indx)) - w + 2
  snp_pos_plot <- snp_pos_total[snp_indx][1:leni]
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
    ]], "[", j = snp_indx)
  }
  # apply rolling window
  df <- lapply(Mres_reduced, function(x) {
    y <- t(apply(x, 1, RcppRoll::roll_sum, n = w))
    return(y)
  })
  df <- do.call(cbind, lapply(lapply(df, function(x) x == 10), colSums,
    na.rm = TRUE
  ))

  df <- data.table::as.data.table(cbind(snp_pos = snp_pos_plot, df))
  df <- data.table::melt(df,
    id.vars = "snp_pos", variable.name = "case_control",
    value.name = "value"
  )
  class(df) <- append("result_snps", class(df))
  return(df)
}

#' print
#'
#' @param x an object
#' @param ... other objects
#'
#' @return print
#' @export
print <- function(x, ...) {
  UseMethod("print")
}

#' print gwas instants
#'
#' @param x object gwas
#'
#' @param ... other objects
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
#'print(snp_data_gds)
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
print.gwas <- function(x, ...) {
  print(paste("object of class", class(x), "with", length(x[[1]]), "samples and", length(x[[2]]), "SNPs"), ...)
}


#' subset an object
#'
#' @param obj object
#'
#' @param ... other variables
#'
#' @export
subset <- function(obj, ...) {
  UseMethod("subset")
}

#' subset gwid object based on snp position
#'
#' @param obj object of class gwid(output of function build_gwid)
#'
#' @param snp_start select starting position of snp, which we want to aggregate.
#' @param snp_end select ending position of snp, which we want to aggregate.
#' @param ... other variables
#'
#' @return the output will be a object(list) of class gwid contains
#' profile object and result_snps object.
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
subset.gwid <- function(obj, snp_start, snp_end, ...) {
  ind_profile <- which(unlist(lapply(obj, inherits, "profile")))
  snp_pos <- obj[["snp_pos"]]
  if (missing(snp_start) & missing(snp_end)) {
    Mres_reduced <- lapply(obj[[ind_profile]], as.matrix)
    class(Mres_reduced) <- "profile"
    output <- list(Mres_reduced = Mres_reduced, snp_pos = snp_pos)
    class(output) <- class(obj)
    return(output)
  }

  if (missing(snp_start)) {
    snp_start <- snp_pos[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_pos[length(snp_pos)]
  }
  snp_indx <- which(snp_pos >= snp_start & snp_pos <= snp_end)
  leni <- diff(range(snp_indx)) + 1
  snp_pos_plot <- snp_pos[snp_indx][1:leni]

  if (length(snp_pos_plot) == length(snp_pos)) {
    Mres_reduced <- lapply(obj[[ind_profile]], as.matrix)
    class(Mres_reduced) <- "profile"

    output <- list(Mres_reduced = Mres_reduced, snp_pos = snp_pos)
    class(output) <- class(obj)
    return(output)
  }
  Mres_reduced <- lapply(lapply(obj[[ind_profile]], "[", j = snp_indx), function(x) {
    x <- as.matrix(x)
    not_zero_index <- which(rowSums(x) != 0)
    y <- x[not_zero_index, ]
  })
  class(Mres_reduced) <- "profile"
  output <- list(Mres_reduced = Mres_reduced, snp_pos = snp_pos_plot)
  class(output) <- class(obj)
  return(output)
}
