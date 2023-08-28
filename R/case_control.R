csv2Rdat <- function(name = "", rep = 3) {
  if (!file.exists(name)) message("File doesn't exist")
  set.seed(1)
  caco <- list()
  output <- list()
  input.data <- utils::read.csv(name)
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
    save(file = paste0(dirname(name),"/",cacos[i], ".Rda"), caco)
    output[[cacos[i]]] <- caco

  }
  output
}


#' Retrieve and load a previously saved file containing a list of case-control data
#'
#'This function loads an RData file that contains an R list containing sample names
#'categorized within different groups. For instance, it includes sample
#'names affiliated with the case group.
#'
#' @param case_control_rda File name in working directory, path to file. A character string giving the name of the case-control
#' file to load.
#'
#' @param ... name of a column (disease name) of csv file.
#'
#' @return The function's outcome will yield an instance of the \strong{caco} class,
#' which takes the form of a list containing names of samples from either the
#' case or control group.
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
case_control <- function(case_control_rda, ...) {
  if (missing(case_control_rda)) {
    stop("provide a case_control list (rda file)")
  }
  substrRight <- function(x, n) {
    substr(x, nchar(x) - n + 1, nchar(x))
  }
  if (substrRight(case_control_rda,3) == "csv") {
    output <- csv2Rdat(case_control_rda)[[...]]
    class(output) <- "caco"
    return(output)
  } else {
    output <- get(load(case_control_rda))
    class(output) <- "caco"
    return(output)
  }
}
