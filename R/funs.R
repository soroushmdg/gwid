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
  output <- list(nas = nas, snp2 = snp2, snps = snps)
  output <- output %>%
    purrr::map(~ dplyr::as_tibble(.x)) %>%
    purrr::map(~ cbind(.x, snp.pos = snp_smp$snp.pos)) %>%
    purrr::map(~ tidyr::pivot_longer(.x, !snp.pos, names_to = "case_control", values_to = "value"))

  class(output) <- append("gwas", class(output))
  return(output)
}


#' @export
plot.gwas <- function(gwas, data, type = "snps") {
  plot_general(gwas, data, type)
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
    # for (k in 1:length(ind[[j]])) {
    #   INDX[ind[[j]][k], 1:2] <- which((snp_smp$smp.id[snp_smp$smp.indx]) %in% IBD[ind[[j]][k], c(1, 3)])
    #   INDX[ind[[j]][k], 3:4] <- which(snp_smp$snp.pos %in% IBD[ind[[j]][k], 6:7])
    #   INDX[ind[[j]][k], 5:6] <- unlist(IBD[ind[[j]][k], c(2, 4)])
    #   LST[[j]][[k]] <- Matrix::Matrix(phased[[INDX[ind[[j]][k], 5]]][INDX[ind[[j]][k], 3]:INDX[ind[[j]][k], 4], INDX[ind[[j]][k], 1]])
    #   mrk1[[j]][[k]] <- as.character(snp_smp$snp.id[INDX[ind[[j]][k], "start"]:INDX[ind[[j]][k], "end"]][which(LST[[j]][[k]][, 1] == 1)])
    # }
  }
  names(ind) <- names(Mres) <- names(LST) <- names(snp_smp$caco)
  output <- list(mrk1 = mrk1, Mres = Mres, LST = LST, INDX = INDX, ind = ind)
  class(output) <- append("raw_gwid", class(output))
  return(output)
}

#' @export
aggregate.raw_gwid <- function(raw_gwid, snp_smp, IBD, ...) {
  res <- res1 <- Subj.id <- list()
  len <- NULL
  IND <- raw_gwid$ind
  LST <- raw_gwid$LST
  Mres <- raw_gwid$Mres
  res <- matrix(0, nr = length(snp_smp$snp.pos), nc = length(snp_smp$caco))
  res1 <- matrix(0, nr = length(snp_smp$snp.pos), nc = length(snp_smp$caco))
  rownames(res1) <- snp_smp$snp.id
  for (j in 1:length(snp_smp$caco)) {
    res[, j] <- apply(Mres[[j]], 2, sum)
    temp <- table(unlist(raw_gwid$mrk1[[j]]))
    res1[names(temp), j] <- temp
  }
  names(IND) <- colnames(res) <- colnames(res1) <- names(snp_smp$caco)
  Subj.id <- (IBD[, c(1, 3)])
  output <- list(LST = LST, Mres = Mres, INDX = raw_gwid$INDX, Subj.id = Subj.id, IND = IND, res = res, res1 = res1)
  output[c(6,7)] <- output[c(6,7)] %>%
    purrr::map(~ dplyr::as_tibble(.x,rownames = NA)) %>%
    purrr::map(~ cbind(.x, snp.pos = snp_smp$snp.pos)) %>%
    purrr::map(~ tidyr::pivot_longer(.x, !snp.pos, names_to = "case_control", values_to = "value"))


  class(output) <- append("gwid", class(output))
  return(output)
}

usethis::use_pipe()



#' @export
plot.gwid <- function(gwid, data, type = "res") {
  plot_general(gwid, data, type)
}


#' @export
count_ibd_window <- function(x, ...){ UseMethod("count_ibd_window") }

#' @export
count_ibd_window.gwid <- function(gwid,snp_smp, w = 10,snp_start, snp_end) {
  if (missing(snp_start)) {
    snp_start <- snp_smp$snp.pos[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_smp$snp.pos[length(snp_smp$snp.pos) - w + 1]
  }


  snp_pos <- snp_smp$snp.pos >= snp_start & snp_smp$snp.pos <= snp_end
  snp_indx <- (which(snp_pos))
  start_end_indx <- range(snp_indx)
  leni <- diff(start_end_indx) - w + 2
  snp_pos_plot <- snp_smp$snp.pos[which(snp_pos)][1:leni]
  if (w >= diff(start_end_indx)) {
    stop("window size should be smaller than number of snps")
  }
  lenj <- length(gwid$Mres)

  Mres <- gwid[["Mres"]]
  Mres <- Mres %>% purrr::map(~as.matrix(.x)) %>% purrr::map(~ dplyr::as_tibble(.x))
  if (w > ncol(Mres[[1]])) stop("window size must be smaller than colums size of the data")
  if (!(is.list(Mres))) stop("data must be of class list!")
  wind_sum <- matrix(0, nrow = leni, ncol = length(Mres))
  for (j in 1:lenj) {
    Mres[[j]] <- Mres[[j]][,snp_indx]
    for (i in 1:leni) {
      wind_sum[i, j] <- sum(apply(Mres[[j]][, i:(w + i - 1), drop = F], 1, sum) == w)
    }
  }
  colnames(wind_sum) <- c("cases", "case1", "case2", "cont1", "cont2", "cont3")
  wind_sum <-wind_sum %>% dplyr::as_tibble() %>% dplyr::bind_cols(snp.pos = snp_smp$snp.pos[snp_indx[1:leni]]) %>%
    tidyr::pivot_longer(!snp.pos, names_to = "case_control", values_to = "value")

  output <- list(window = wind_sum)

  class(output) <- append("count_ibd", class(output))
  return(output)
}

#' @export
plot.count_ibd <- function(count_ibd, snp_smp, type = "window") {
  plot_general( statistics = count_ibd,data = snp_smp , type)
}


# #' @export
# haplo_str_window <- function(x, ...){UseMethod("haplo_str_window") }

#' @export
haplo_str_window <- function(gwid, snp_smp, phased, w = 10, snp_start, snp_end) {
  if (missing(snp_start)) {
    snp_start <- snp_smp$snp.pos[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_smp$snp.pos[length(snp_smp$snp.pos) - w + 1]
  }
  snp_pos <- snp_smp$snp.pos >= snp_start & snp_smp$snp.pos <= snp_end
  snp_indx <- (which(snp_pos))
  start_end_indx <- range(snp_indx)
  leni <- diff(start_end_indx) - w + 2
  snp_pos_plot <- snp_smp$snp.pos[which(snp_pos)][1:leni]
  if (w >= diff(start_end_indx)) {
    stop("window size should be smaller than number of snps")
  }
  lenj <- length(gwid$Mres)
  structures <- vector(mode = "list", length = lenj)
  names(structures) <- names(gwid$Mres)

  structures_gwas <- vector(mode = "list", length = lenj)
  names(structures_gwas) <- names(snp_smp$caco)

  for (j in 1:lenj) {
    caco_j <- which(colnames(phased[[1]]) %in% snp_smp$caco[[j]])
    structures[[j]] <- vector(mode = "list", length = leni)
    structures_gwas[[j]] <- vector(mode = "list", length = leni)
    for (i in 1:leni) {
      structures_gwas[[j]][[i]] <- as.vector(apply(cbind(
        phased[[1]][snp_indx[i]:(w + snp_indx[i] - 1),
          caco_j,
          drop = F
        ],
        phased[[2]][snp_indx[i]:(w + snp_indx[i] - 1),
          caco_j,
          drop = F
        ]
      ), 2, paste, collapse = ""))

      ibd_reg_count <- which(apply(gwid$Mres[[j]][, snp_indx[i]:(w + snp_indx[i] - 1), drop = F], 1, sum) == w)
      for (k in 1:length(ibd_reg_count)) {
        indk <- gwid$IND[[j]][ibd_reg_count[k]]
        structures[[j]][[i]][k] <- paste0(gwid$LST[[j]][[ibd_reg_count[k]]][((snp_indx[i] - gwid$INDX[indk, "start"]) + (1:w))], collapse = "")
      }
    }
  }
  names(structures) <- c("cases", "case1", "case2", "cont1", "cont2", "cont3")
  output <- list(structures = structures, structures_gwas = structures_gwas, snp_pos_plot = snp_pos_plot)
  class(output) <- append("haplo_str", class(output))
  return(output)
}

#' @export
haplo_str_window2 <- function(gwid, ibd_data, snp_smp, phased, w = 10, snp_start = 2.5 * 10^7, snp_end = 2.6 * 10^7) {
  #browser()
  if (missing(snp_start)) {
    snp_start <- snp_smp$snp.pos[1]
  }
  if (missing(snp_end)) {
    snp_end <- snp_smp$snp.pos[length(snp_smp$snp.pos)]
  }
  snp_pos <- snp_smp$snp.pos >= snp_start & snp_smp$snp.pos <= snp_end
  snp_indx <- which(snp_pos)
  start_end_indx <- range(snp_indx)
  leni <- diff(start_end_indx) - w + 2
  snp_pos_plot <- snp_smp$snp.pos[which(snp_pos)][1:leni]
  if (w >= diff(start_end_indx)) {
    stop("window size should be smaller than number of snps")
  }

  # apply rolling window
  df <- gwid$Mres %>%
    purrr::map(~ as.matrix(.x)) %>%
    purrr::map(~ dplyr::as_tibble(.x)) %>%
    purrr::map(~ dplyr::select(.x, snp_indx)) %>%
    purrr::map(~ (t(apply(.x, 1, RcppRoll::roll_sum, w)))) %>%
    purrr::map(~ dplyr::as_tibble(.x))

  # adding haplotype and sample information
  for (i in 1:length(df)) {
    haplotype <- ibd_data[gwid$IND[[i]], 2] %>%
      dplyr::as_tibble() %>%
      dplyr::rename("haplo_index" = V2)
    smp <- ibd_data[gwid$IND[[i]], 1] %>%
      dplyr::as_tibble() %>%
      dplyr::rename("smp" = V1)
    df[[i]] <- dplyr::bind_cols(haplotype, smp, df[[i]])
  }

  # identify sample locations in phased data
  sample_indexing <- dplyr::tibble(smp = colnames(phased[[1]])) %>% tibble::rowid_to_column(var = "smp_index")


  # locate in-ibd regions
  df <- df %>%
    purrr::map(~ tidyr::pivot_longer(.x, !c(haplo_index, smp), values_to = "sum_window")) %>%
    purrr::map(~ dplyr::filter(.x, sum_window == w)) %>%
    purrr::map(~ dplyr::mutate(.x, window_number = as.numeric(stringr::str_sub(name, 2, -1)))) %>%
    purrr::map(~ dplyr::select(.x, haplo_index, smp, window_number)) %>%
    dplyr::bind_rows(.id = "case_control") %>%
    dplyr::left_join(x = ., y = sample_indexing, by = "smp") %>%
    dplyr::select(case_control, haplo_index, smp, smp_index, window_number) %>%
    dplyr::arrange(window_number)


  result <- Matrix(NA, nr = nrow(df), nc = w)
  my_phase1 <- phased[[1]][start_end_indx[1]:nrow(phased[[1]]), ]
  my_phase2 <- phased[[2]][start_end_indx[1]:nrow(phased[[2]]), ]

  for (i in 1:w) {
    result[df$haplo_index == 1, i] <- my_phase1[as.matrix(df[df$haplo_index == 1, c(5, 4)])]
    result[df$haplo_index == 2, i] <- my_phase2[as.matrix(df[df$haplo_index == 2, c(5, 4)])]
    df[, 5] <- df[, 5] + 1
  }

  df[, 5] <- df[, 5] - w
  res <- dplyr::bind_cols(df, structures = apply(result, 1, paste, collapse = ""))
  output <- list(haplo_str = res,snp_pos = snp_pos_plot)
  return(output)
}

#' @export
haplo_freq2 <- function(hap_str){
  df <-  hap_str[[1]] %>%
    dplyr::group_by(window_number,case_control,structures) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::mutate(window_number= factor(window_number, levels = unique(window_number)))
  # haplotype frequency
  df <- split(df,f = df$window_number) %>%
    purrr::map(~ dplyr::mutate(.x,structures = factor(structures,levels = unique(structures)) )) %>%
    purrr::map(~ tidyr::complete(data = .x, structures, fill = list(n = 0))) %>%
    dplyr::bind_rows()

  # haplotype most frequency
  df1 <- df %>%
    dplyr::group_by(window_number,case_control) %>%
    dplyr::summarise(most_freq = max(n)) %>%
    dplyr::group_by(case_control) %>%
    dplyr::group_split() %>%
    purrr::map(~ dplyr::bind_cols(.x,snp_pos = hap_str2[[2]] )) %>%
    dplyr::bind_rows()

  # number of ibd in each window
  df2 <- df %>%
    dplyr::group_by(window_number,case_control) %>%
    dplyr::summarise(value = sum(n)) %>%
    dplyr::group_by(case_control) %>%
    dplyr::group_split() %>%
    purrr::map(~ dplyr::bind_cols(.x,snp_pos = hap_str2[[2]] )) %>%
    dplyr::bind_rows()

  output <- list(hap_freq = df, hap_most_freq = df1, count_ibd = df2)
  output <- output %>% purrr::map(~ dplyr::mutate(.x,window_number = as.numeric(window_number)))
  return(output)

}

#' @export
haplo_freq_window <- function(x, ...){UseMethod("haplo_freq_window") }

#' @export
haplo_freq_window.haplo_str <- function(hap_str) {
  output <- list()
  for (i in 1:2) {
    freq <- matrix(unlist(lapply(lapply(unlist(hap_str[[i]], recursive = F), table), max)), ncol = length(hap_str[[i]]), byrow = F)
    colnames(freq) <- c("cases", "case1", "case2", "cont1", "cont2", "cont3")
    mywindow <- rep(1:nrow(freq))
    snp.pos <- hap_str$snp_pos_plot
    freq <- cbind(mywindow, snp.pos, freq)
    freq <- dplyr::as_tibble(freq)
    freq <- freq %>% tidyr::pivot_longer(!c(mywindow, snp.pos), names_to = "case_control", values_to = "structure")
    # unlist haplo structures
    new.l <- rapply(hap_str[[i]], function(x) paste(x, collapse = "|"), how = "replace")
    dt <- data.table::rbindlist(new.l)
    dt <- data.table::transpose(dt)
    colnames(dt) <- names(hap_str[[i]])
    dt <- cbind(mywindow, snp.pos, dt)
    col_names <- colnames(dt)
    dt <- dt %>% tidyr::pivot_longer(!c(mywindow, snp.pos), names_to = "case_control", values_to = "structure")
    dt <- dt %>% tidyr::separate_rows(structure, sep = "\\|")
    hap_freq <- dt %>%
      dplyr::group_by(mywindow, case_control, structure) %>%
      dplyr::summarise(n = dplyr::n())
    output[[i]] <- list(hap_most_freq = freq, hap_freq = hap_freq)
  }
  names(output) <- c("gwid", "gwas")
  class(output) <- append("haplo_freq", class(output))
  return(output)
}



#' @export
plot_general <- function(df, type = "aggregate_gwid", nwin) {
  if (type %in% c("aggregate_gwas", "aggregate_gwid", "count_ibd_window","haplo_most_freq")) {
    p <- df %>% ggplot2::ggplot(ggplot2::aes(x = snp.pos, y = value)) +
      ggplot2::geom_line(ggplot2::aes(color = case_control), size = .6) +
      ggplot2::scale_x_continuous("snp position", labels = scales::label_number_si()) +
      ggplot2::scale_y_continuous("sum of type") +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 10))

    fig <- plotly::ggplotly(p)
  }

  if (type == "structure") {
      p <- df %>%
        dplyr::mutate(
          mywindow = factor(mywindow, levels = unique(mywindow)),
          structure = factor(structure, levels = unique(structure))
        ) %>%
        dplyr::filter(mywindow == nwin) %>%
        ggplot2::ggplot(ggplot2::aes(x = structure, y = n)) +
        ggplot2::geom_line(ggplot2::aes(color = case_control, group = case_control), size = .6) +
        #ggplot2::scale_x_discrete("haplotype structures") +
        ggplot2::scale_y_discrete("sum of structures") +
        ggplot2::theme(legend.title = ggplot2::element_text(size = 10)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggplot2::facet_wrap(ggplot2::vars(mywindow), scales = "free_x")
      fig <- plotly::ggplotly(p)
    }


  if (type == "structure_v2") {
    fig <- list()

      p <- df %>%
        dplyr::mutate(
          mywindow = factor(mywindow, levels = unique(mywindow)),
          structure = factor(structure, levels = unique(structure))
        ) %>%
        dplyr::filter(mywindow == nwin) %>%
        ggplot2::ggplot(ggplot2::aes(x = case_control, y = n)) +
        ggplot2::geom_line(ggplot2::aes(color = structure, group = structure), size = .6) +
        #ggplot2::scale_x_discrete("haplotype structures") +
        ggplot2::scale_y_discrete("sum of structures") +
        ggplot2::theme(legend.title = ggplot2::element_text(size = 10)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
      fig <- plotly::ggplotly(p)
    }

  fig
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

