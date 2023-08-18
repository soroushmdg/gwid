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
#' @param ... other variables
#'
#' @return an interactive line plot of result_snps for each case control subjects.
#'
#' @export
plot.result_snps <- function(x, y = NA, title, snp_start, snp_end, ly = TRUE, ...) {
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
    ggplot2::geom_line(ggplot2::aes(color = case_control), size = .6) +
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
#' @param ... other variables
#'
#' @return an interactive line plot of test_snps objects for each case control subjects.
#'
#' @export
plot.test_snps <- function(x, y = NA, title, snp_start, snp_end, ly =TRUE, ...) {
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
    ggplot2::geom_line(ggplot2::aes(color = case_control), size = .6) +
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
#' @param ... other variables
#' @export
plot.haplotype_structure_frequency <- function(x, y = NA, type = c("version1", "version2"), nwin, ly = TRUE, ...) {
  type <- match.arg(type)
  structures <- window_number <- n <- NULL
  if (any(!is.na(y))){
    x <- x[case_control %in% y]
  }
  if (type == "version1") {
    p <- ggplot2::ggplot(x[window_number == nwin, ], ggplot2::aes(x = structures, y = n, frame = window_number)) +
      ggplot2::geom_line(ggplot2::aes(color = case_control, group = case_control), size = .6) +
      ggplot2::scale_x_discrete("haplotypes", label = function(x) strtrim(x, 12)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggplot2::labs(
        title = "haplotype frequency version 1", y = "number of structures",
        fill = "case_control"
      )
  } else {
    p <- ggplot2::ggplot(x[window_number == nwin, ], ggplot2::aes(x = case_control, y = n, frame = window_number)) +
      ggplot2::geom_line(ggplot2::aes(color = structures, group = structures), size = .6) +
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
#' @param title title of the plot.
#' @param ... optional argument of \code{plot}
#'
#' @export
plot.haplotype_frequency <- function(x, y = NA, plot_type = c("haplotype_structure_frequency", "result_snps"), type = c("version1", "version2"), nwin, title, ...) {
  plot_type <- match.arg(plot_type)
  if (plot_type == "haplotype_structure_frequency") {
    type <- match.arg(type)
    plot(x[[which(unlist(lapply(x, inherits, "haplotype_structure_frequency")))]], y, type = type, nwin = nwin)
  } else {
    if (!missing(title)) {
      plot(x[[which(unlist(lapply(x, inherits, "result_snps")))]], y, title = title, ...)
    } else {
      plot(x[[which(unlist(lapply(x, inherits, "result_snps")))]], y, ...)
    }
  }
}
