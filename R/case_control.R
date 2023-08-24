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
