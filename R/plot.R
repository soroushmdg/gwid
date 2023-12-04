#'  Line plot of result_snps objects
#'
#' @param x An object of class result_snps.
#' @param y default value is NA, if specified it should be a vector of names of
#' subject groups i.e. y = c("case","control")
#'
#' @param title title of the plot.
#' @param snp_start select starting position of snps.
#' @param snp_end select ending position of snps.
#' @param ly if TRUE, we have a plotly object and if it is false plot is going to be
#' a ggplot object.
#' @param line_size geom_line size
#' @param ... other variables
#'
#' @return an interactive line plot of result_snps for each case control subjects.
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
plot.result_snps <- function(x, y = NA, title, snp_start, snp_end, ly = TRUE, line_size = .6, ...) {
  snp_pos <- value <- NULL
  if (!missing(snp_start) & !missing(snp_end)) {
    x <- x[snp_pos >= snp_start & snp_pos <= snp_end]
  }
  if (!missing(snp_start) & missing(snp_end)) {
    x <- x[snp_pos >= snp_start]
  }

  if (missing(snp_start) & !missing(snp_end)) {
    x <- x[snp_pos <= snp_end]
  }
  if (any(!is.na(y))){
    x <- x[case_control %in% y]
  }
  p <- ggplot2::ggplot(x, ggplot2::aes(x = snp_pos, y = value)) +
    ggplot2::geom_line(ggplot2::aes(color = case_control), size = line_size) +
    ggplot2::scale_x_continuous("snp position",
      labels = paste0(round(quantile(x$snp_pos, seq(0, 1, length.out = 5)) / 10^6), "M"),
      breaks = quantile(x$snp_pos, seq(0, 1, length.out = 5))
    )

  if (missing(title)) {
    p <- p + ggplot2::labs(
      y = "count",
      fill = "case_control"
    )
  } else {
    p <- p + ggplot2::labs(
      title = title, y = "count",
      fill = "case_control"
    )
  }
  p <- p + ggplot2::theme(legend.title = ggplot2::element_text(size = 10)) + ggplot2::theme(
    axis.ticks = ggplot2::element_line(colour = "black"),
    panel.grid.major = ggplot2::element_line(
      colour = "antiquewhite",
      linetype = "dashed"
    ), panel.grid.minor = ggplot2::element_line(
      colour = "antiquewhite",
      linetype = "dashed"
    ), axis.text = ggplot2::element_text(colour = "black"),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.title = ggplot2::element_text(size = 12),
    panel.background = ggplot2::element_rect(
      fill = NA,
      colour = NA
    ),
    plot.background = ggplot2::element_rect(colour = NA),
    legend.key = ggplot2::element_rect(fill = NA),
    legend.background = ggplot2::element_rect(fill = NA)
  ) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(colour = NA),
      plot.background = ggplot2::element_rect(colour = NA)
    ) +
    ggplot2::scale_colour_brewer(type = "qual", palette = 2)
  if (ly == TRUE){
    p <- plotly::ggplotly(p, ...)
    return(p)
  } else {
    return(p)
  }
}

#'  Line plot of gwid objects
#'
#' @param x An object of class gwid. Output of \code{build_gwid} function.
#' @param y default value is NA, if specified it should be a vector of names of
#' subject groups i.e. y = c("case","control")
#'
#' @param title title of the plot.
#' @param plot_type either \dQuote{result_snps} or \dQuote{profile}.
#' @param reference reference group of subjects in which we want to have profile plot.
#' @param ... if plot_type is \dQuote{result_snps} it is optional argument of \code{plot}.
#' if plot_type is \dQuote{profile} we can subset plot based on snp_start and snp_end locations.
#'
#' @return if plot_type is \dQuote{result_snps} an interactive line plot of result_snps for each case control subjects.
#' if plot_type is \dQuote{profile} an interactive profile plot of identity by descent subjects in subset of locations.
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
plot.gwid <- function(x, y = NA, title = "number of IBD in each snp", plot_type = c("result_snps", "profile"), reference, ...) {
  plot_type <- match.arg(plot_type)
  if (plot_type == "profile") {
    dots <- list(...)
    invisible(ifelse("snp_start" %in% names(dots), assign("snp_start", dots[["snp_start"]]), assign("snp_start", 0)))
    invisible(ifelse("snp_end" %in% names(dots), assign("snp_end", dots[["snp_end"]]), assign("snp_end", Inf)))
    snp_pos_ind <- which(x$snp_pos >= snp_start & x$snp_pos <= snp_end)
    z <- as.matrix(x[[which(unlist(lapply(x, inherits, "profile")))]][[reference]])[, snp_pos_ind]
    z <- z[which(rowSums(z) != 0), ]
    snp_pos <- x$snp_pos[snp_pos_ind]

    plotly::plot_ly(x = snp_pos, z = z, type = "heatmap", colors = "Greys", showscale = FALSE) |>
      plotly::layout(title = reference, xaxis = list(title = "Positions"), yaxis = list(title = "Profiles"))
  } else {
    if (sum(unlist(lapply(x, inherits, "result_snps"))) == 0) {
      x$res <- extract(x)
    }

    plot(x[[which(unlist(lapply(x, inherits, "result_snps")))]], y,
      title = title, ...
    )
  }
}

#'  Line plot of gwas objects
#'
#' @param x object of class gwas.
#' @param y default value is NA, if specified it should be a vector of names of
#' subject groups i.e. y = c("case","control")
#' @param title title of the plot.
#' @param ... optional argument of \code{plot}
#'
#' @return an interactive line plot of gwas objects for each case control subjects.
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
plot.gwas <- function(x, y = NA, title = "number of snps", ...) {
  if (sum(unlist(lapply(x, inherits, "result_snps"))) == 0) {
    x$snps <- extract(x)
  }

  plot(x[[which(unlist(lapply(x, inherits, "result_snps")))]], y,
    title = title, ...
  )
}



#' Line plot of test_snps objects
#'
#' @param x an object of class test_snps.
#' @param y default value is NA, if specified it should be a vector of names of
#' subject groups i.e. y = c("case","control")
#' @param title title of the plot.
#' @param snp_start select starting position of snps.
#' @param snp_end select ending position of snps.
#' @param ly if `TRUE`, we have a `plotly` object and if it is `FALSE` plot is going to be
#' a `ggplot` object.
#' @param line_size geom_line size
#' @param ... other variables
#'
#' @return an interactive line plot of test_snps objects for each case control subjects.
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
plot.test_snps <- function(x, y = NA, title, snp_start, snp_end, ly =TRUE, line_Size = .6, ...) {
  snp_pos <- value <- NULL
  if (!missing(snp_start) & !missing(snp_end)) {
    x <- x[snp_pos >= snp_start & snp_pos <= snp_end]
  }
  if (!missing(snp_start) & missing(snp_end)) {
    x <- x[snp_pos >= snp_start]
  }

  if (missing(snp_start) & !missing(snp_end)) {
    x <- x[snp_pos <= snp_end]
  }
  if (any(!is.na(y))){
    x <- x[case_control %in% y]
  }
  p <- ggplot2::ggplot(x, ggplot2::aes(x = snp_pos, y = value)) +
    ggplot2::geom_line(ggplot2::aes(color = case_control), size = line_Size) +
    ggplot2::scale_x_continuous("snp position",
      labels = paste0(round(quantile(x$snp_pos, seq(0, 1, length.out = 5)) / 10^6), "M"),
      breaks = quantile(x$snp_pos, seq(0, 1, length.out = 5))
    ) +
    ggplot2::scale_y_continuous(trans = "log10")
  if (!missing(title)) {
    p <- p + ggplot2::labs(
      title = title, y = "p_values",
      fill = "case_control"
    )
  } else {
    p <- p + ggplot2::labs(
      y = "p_values",
      fill = "case_control"
    )
  }
  p <- p + ggplot2::theme(legend.title = ggplot2::element_text(size = 10)) + ggplot2::theme(
    axis.ticks = ggplot2::element_line(colour = "black"),
    panel.grid.major = ggplot2::element_line(
      colour = "antiquewhite",
      linetype = "dashed"
    ), panel.grid.minor = ggplot2::element_line(
      colour = "antiquewhite",
      linetype = "dashed"
    ), axis.text = ggplot2::element_text(colour = "black"),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.title = ggplot2::element_text(size = 12),
    panel.background = ggplot2::element_rect(
      fill = NA,
      colour = NA
    ),
    plot.background = ggplot2::element_rect(colour = NA),
    legend.key = ggplot2::element_rect(fill = NA),
    legend.background = ggplot2::element_rect(fill = NA)
  ) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(colour = NA),
      plot.background = ggplot2::element_rect(colour = NA)
    ) +
    ggplot2::scale_colour_brewer(type = "qual",  palette = 2)
  if (ly == TRUE){
    p <- plotly::ggplotly(p, ...)
    return(p)
  } else {
    return(p)
  }


}


#' Two type of line plots for haplotype_structure_frequency objects .
#'
#' @param x an object of class haplotype_structure_frequency
#' @param y default value is NA, if specified it should be a vector of names of
#' subject groups i.e. y = c("case","control")
#' @param type either \dQuote{version1} or \dQuote{version2}
#' @param nwin window number
#' @param ly if `TRUE`, we have a `plotly` object and if it is `FALSE` plot is going to be
#' a `ggplot` object.
#' @param line_size geom_line size
#' @param ... other variables
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
plot.haplotype_structure_frequency <- function(x, y = NA, type = c("version1", "version2"), nwin, ly = TRUE, line_size=.6,...) {
  type <- match.arg(type)
  structures <- window_number <- n <- NULL
  if (any(!is.na(y))){
    x <- x[case_control %in% y]
  }
  if (type == "version1") {
    p <- ggplot2::ggplot(x[window_number == nwin, ], ggplot2::aes(x = structures, y = n, frame = window_number)) +
      ggplot2::geom_line(ggplot2::aes(color = case_control, group = case_control), size = line_size) +
      ggplot2::scale_x_discrete("haplotypes", label = function(x) strtrim(x, 12)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggplot2::labs(
        title = "haplotype frequency version 1", y = "number of structures",
        fill = "case_control"
      )
  } else {
    p <- ggplot2::ggplot(x[window_number == nwin, ], ggplot2::aes(x = case_control, y = n, frame = window_number)) +
      ggplot2::geom_line(ggplot2::aes(color = structures, group = structures), size = line_size) +
      ggplot2::scale_x_discrete("structures") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggplot2::labs(
        title = "haplotype frequency version 2", y = "number of structures",
        fill = "case_control"
      )
  }
  p <- p + ggplot2::theme(legend.title = ggplot2::element_text(size = 10)) + ggplot2::theme(
    axis.ticks = ggplot2::element_line(colour = "black"),
    panel.grid.major = ggplot2::element_line(
      colour = "antiquewhite",
      linetype = "dashed"
    ), panel.grid.minor = ggplot2::element_line(
      colour = "antiquewhite",
      linetype = "dashed"
    ), axis.text = ggplot2::element_text(colour = "black"),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.title = ggplot2::element_text(size = 12),
    panel.background = ggplot2::element_rect(
      fill = NA,
      colour = NA
    ),
    plot.background = ggplot2::element_rect(colour = NA),
    legend.key = ggplot2::element_rect(fill = NA),
    legend.background = ggplot2::element_rect(fill = NA)
  ) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(colour = NA),
      plot.background = ggplot2::element_rect(colour = NA)
    )
  if (type != "version2") {
    p <- p + ggplot2::scale_colour_brewer(type = "qual", palette = 2)
  }
  if (ly == TRUE){
    p <- plotly::ggplotly(p)
    return(p)
  } else {
    return(p)
  }

}


#' Line plot of haplotype_frequency object
#'
#' @param x an object of class haplotype_frequency
#' @param y default value is `NA`, if specified it should be a vector of names of
#' subject groups i.e. `y = c("case","control")`
#' @param plot_type either \dQuote{result_snps} or \dQuote{"haplotype_structure_frequency"}
#' @param type either \dQuote{version1} or \dQuote{version2} when plot_type is \dQuote{"haplotype_structure_frequency"}
#' @param nwin window number
#' @param ly if TRUE, we have a plotly object and if it is false plot is going to be
#' a ggplot object.
#' @param title title of the plot.
#' @param line_size geom_line size
#' @param ... optional argument of \code{plot}
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
plot.haplotype_frequency <- function(x, y = NA, plot_type = c("haplotype_structure_frequency", "result_snps"), type = c("version1", "version2"),ly=TRUE, nwin, title,line_size = .6, ...) {
  plot_type <- match.arg(plot_type)
  if (plot_type == "haplotype_structure_frequency") {
    type <- match.arg(type)
    plot(x[[which(unlist(lapply(x, inherits, "haplotype_structure_frequency")))]], y, type = type, nwin = nwin,ly = ly,line_size=line_size)
  } else {
    if (!missing(title)) {
      plot(x[[which(unlist(lapply(x, inherits, "result_snps")))]], y, title = title, ly=ly, ...)
    } else {
      plot(x[[which(unlist(lapply(x, inherits, "result_snps")))]], y, ly = ly, ...)
    }
  }
}
