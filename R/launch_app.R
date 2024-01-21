shiny_ui <-
  fluidPage(
    tags$head(tags$style(HTML("body { max-width: 1250px !important; }"))),
    titlePanel("gwid"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        tags$head(tags$style(type = "text/css", ".well { max-width: 300px; }")),
        uiOutput("f_choice", width = "400px"),
        sliderInput("chr", "Select Chromosome",
          min = 1, max = 22,
          value = 22, step = 1, width = "400px"
        ),
        uiOutput("sh_plot", width = "400px"),
        uiOutput("res_choice", width = "300px"),
        conditionalPanel(
          condition = "input.sh_plot=='rohm' || input.sh_plot=='rohs'",
          sliderInput("wroh", "ROH Win",
            min = 2, max = 20, step = 1, value = 5,
            width = "200px"
          )
        ),
        uiOutput("download.b"),
        tags$hr(style = "border-color: red;", width = "150px"),
        conditionalPanel(
          condition = "input.sh_plot=='rohm' || input.sh_plot=='rohs'",
          checkboxInput("roh_test", "ROH (McNemar.Test)",
            value = FALSE,
            width = "400px"
          ),
          column(8, uiOutput("mcnmr_choice", width = "240px")),
          column(4, uiOutput("wmcn", width = "200px")),
          checkboxInput("roh_perm", "ROH (Permutation.Test)",
            value = FALSE, width = "400px"
          ),
          column(6, uiOutput("nperm", width = "200px")),
          column(6, uiOutput("runtestr", width = "200px"))
        ),
        conditionalPanel(
          condition = " input.Panel=='IBD' ",
          checkboxInput("genome_ibd", "Genome (vs IBD reg.)",
            value = FALSE, width = "400px"
          )
        ),
        # uiOutput("skip.rep", width = "400px"),
        tags$hr(style = "border-color: red;", width = "150px"),
        conditionalPanel(
          condition = "input.sh_plot != 'rohs' && input.sh_plot != 'rohm'",
          checkboxInput("fisher", "Fisher Exact (G-) Test",
            value = FALSE,
            width = "400px"
          ),
          conditionalPanel(
            condition = "input.fisher",
            column(8, uiOutput("fisher_choice")),
            column(4, uiOutput("fisher_side"))
          ),
          checkboxInput("perm_gap", "permutation test",
            value = FALSE,
            width = "400px"
          ),
          conditionalPanel(
            condition = "input.perm_gap",
            column(8, uiOutput("nperm_gwid"))
          ),
          conditionalPanel(
            condition = "input.perm_gap || input.fisher",
            uiOutput("runtest", width = "400px")
          ),
          conditionalPanel(
            condition = "input.Panel=='IBD' && (input.fisher || input.perm_gap) ",
            uiOutput("runtestw", width = "400px")
          )
        ),
      ),
      mainPanel(
        width = 9,
        tabsetPanel(
          id = "Panel", type = "tabs",
          tabPanel("Select Region",
            value = "nIBD",
            column(12, uiOutput("snp_selected", align = "center"),
              style = "color:red;"
            ),
            column(12, uiOutput("pos", width = "900px")),
            column(12, plotOutput("reg_out", height = 90, width = 900)),
            column(12, plotly::plotlyOutput("res_out",
              height = 300,
              width = 900
            )),
            column(12, plotly::plotlyOutput("test_out",
              height = 300,
              width = 900
            ))
          ),
          tabPanel("IBD (Window Based)",
            value = "IBD",
            column(12, uiOutput("snp_selected2", align = "center"),
              style = "color:red;"
            ),
            plotOutput("reg_out2", height = 100, width = 900),
            column(3, uiOutput("s_pos", width = "200px")),
            column(3, uiOutput("e_pos", width = "200px")),
            column(3, uiOutput("runit", width = "200px")),
            column(3, uiOutput("wsize", width = "200px")),
            column(12, plotly::plotlyOutput("IBD_out",
              height = 300,
              width = 900
            )),
            column(12, uiOutput("snp_win", width = "900px")),
            column(12, plotly::plotlyOutput("test_outw",
              height = 300,
              width = 900
            ))
          )
        )
      )
    )
  )


shiny_server <- function(input, output, session) {
  ival <- reactiveVal(c(0, 0))

  values <- reactiveValues()
  data_path <- reactiveVal()


  # changing data type, event -----------------------------------------------


  observeEvent(input$f_choice,
    {
      if (is.null(input$f_choice)) {
        return()
      }
      data_path(paste0(data_folder_address, "/"))
      values$caco <- case_control(case_control_rda = paste0(data_path(), input$f_choice))
      values$pieces <- build_gwas(
        gds_data = paste0(data_path(), "chr", input$chr, ".gds"),
        caco = values$caco, gwas_generator = TRUE)
      rng <- range(values$pieces[["snp.pos"]])
      updateSliderInput(session, "pos", min = rng[1], max = rng[2], value = rng)
      values$ag.pieces <- values$pieces$snps
      values$myphase <- build_phase(
        phased_vcf = paste0(data_path(), "chr", input$chr, ".vcf"),
        caco = values$caco)
      values$myregion2 <- build_gwid(
        ibd_data = paste0(data_path(), "chr", input$chr, ".ibd"),
        gwas = values$pieces, gwid_generator = TRUE)
      values$ibd <- values$myregion2$ibd
      updateSelectInput(session, "res_choice",
                        choices = names(values$caco),
                        selected = names(values$caco)[1])
    },
    priority = 2
  )


  # changing chromosome, event ----------------------------------------------


  observeEvent(input$chr,
    {
      if (is.null(input$f_choice)) {
        return()
      }

      if (!file.exists(paste0(data_path(), "chr", input$chr, ".gds"))) {
        showNotification("There is no data associated to this chromosome",
          type = "message"
        )
        print("There is no data associated to this chromosome")
        return()
      }
      data_path(paste0(data_folder_address, "/"))
      values$caco <- case_control(case_control_rda = paste0(data_path(),
                                                            input$f_choice))
      values$pieces <- build_gwas(
        gds_data = paste0(
          data_path(),
          "chr", input$chr, ".gds"
        ),
        caco = values$caco, gwas_generator = TRUE
      )
      rng <- range(values$pieces[["snp.pos"]])

      if (ival()[1] < rng[1]) {
        ival(c(rng[1], ival()[2]))
      }
      if (ival()[2] > rng[2]) {
        ival(c(ival()[1], rng[2]))
      }



      updateSliderInput(session, "pos",
        min = rng[1], max = rng[2],
        value = ival()
      )
      values$ag.pieces <- values$pieces$snps
      values$myphase <- build_phase(
        phased_vcf = paste0(
          data_path(), "chr",
          input$chr, ".vcf"
        ),
        caco = values$caco
      )
      values$myregion2 <- build_gwid(
        ibd_data = paste0(
          data_path(), "chr",
          input$chr, ".ibd"
        ),
        gwas = values$pieces, gwid_generator = TRUE
      )
      values$ibd <- values$myregion2$ibd
    },
    priority = 1
  )

  observeEvent(input$roh_test, {
    if (input$roh_test) {
      updateCheckboxInput(session, "roh_perm", value = FALSE)
      updateCheckboxInput(session, "fisher", value = FALSE)
      updateCheckboxInput(session, "perm_gap", value = FALSE)
    }
  })

  observeEvent(input$roh_perm, {
    if (input$roh_perm) {
      updateCheckboxInput(session, "roh_test", value = FALSE)
      updateCheckboxInput(session, "fisher", value = FALSE)
      updateCheckboxInput(session, "perm_gap", value = FALSE)
    }
  })

  observeEvent(input$fisher, {
    if (input$fisher) {
      updateCheckboxInput(session, "roh_test", value = FALSE)
      updateCheckboxInput(session, "roh_perm", value = FALSE)
      updateCheckboxInput(session, "perm_gap", value = FALSE)
    }
  })

  observeEvent(input$perm_gap, {
    if (input$perm_gap) {
      updateCheckboxInput(session, "roh_test", value = FALSE)
      updateCheckboxInput(session, "roh_perm", value = FALSE)
      updateCheckboxInput(session, "fisher", value = FALSE)
    }
  })

  # RA-withmap, RA, ... -----------------------------------------------------


  output$f_choice <- renderUI({
    f_choices <- list.files(
      path = data_folder_address, pattern = "\\.rda$",
      ignore.case = TRUE
    )
    selectInput("f_choice", "Select a Data",
      choices = f_choices,
      selected = f_choices[1], width = "400px"
    )
  })

  # plot types --------------------------------------------------------------


  output$sh_plot <- renderUI({
    if (input$Panel != "IBD") {
      p.choices <- list(
        "Genotype" = c(
          "SNPs in IBD regions" = "snp_ibd",
          "Profile Plot" = "profile"
        ),
        "Basic GWAS" = c("Number of SNPs" = "snp_cnt")
      )
      if (!is.null(input$f_choice)) {
        p.choices <- append(
          p.choices,
          list(
            ROH =
              c(
                "ROH Mean" =
                  "rohm",
                "ROH Sum" =
                  "rohs"
              )
          ), 1
        )
      }
    } else {
      p.choices <- list(
        "Genotype" = c("SNPs in IBD window" = "snp_win"),
        "Phased" = c(
          "# Most Freq. Hapl." = "hap_win",
          "Hapl. v1" = "hap_1", "Hapl. v2" = "hap_2"
        )
      )
      names(p.choices[[2]]) <- paste(names(p.choices[[2]]), "in IBD win.")
    }
    selectInput("sh_plot", "Select a Plot",
      choices = p.choices,
      selected = ifelse(input$Panel != "IBD", "snp_ibd", "snp_win"),
      width = "400px"
    )
  })

  # visualization 1 ---------------------------------------------------------
  output$res_out <- plotly::renderPlotly({
    if (is.null(values$myregion2)) {
      return()
    }
    if (input$Panel != "nIBD") {
      return()
    }

    if (input$sh_plot %in% c("snp_win", "hap_win", "hap_1", "hap_2")) {
      return()
    }
    if (!setequal(ival(), c(0, 0))) {
      snp_pos <- values$pieces$snp.pos >= ival()[1] &
        values$pieces$snp.pos <= ival()[2]


      if (input$sh_plot == "snp_ibd") {
        p <- plot(values$myregion2,
          title = "gwid", snp_start = ival()[1],
          snp_end = ival()[2]
        )
      } else if (input$sh_plot == "profile") {
        if (is.null(input$res_choice)) {
          return()
        }
        len <- sum(values$pieces[["snp.pos"]] >= input$pos[1] &
          values$pieces[["snp.pos"]] <= input$pos[2])

        if (len >= 5000) {
          return()
        }
        values$myprofile <- subset(values$myregion2,
          snp_start = ival()[1],
          snp_end = ival()[2]
        )
        p <- plot(values$myprofile,
          plot_type = "profile",
          reference = input$res_choice
        )
      } else if (input$sh_plot == "snp_cnt") {
        ag.pieces <- values$ag.pieces
        ag.pieces <- ag.pieces[snp_pos >= ival()[1] & snp_pos <= ival()[2], ]
        p <- plot(ag.pieces, title = "gwas")
      } else if (input$sh_plot == "rohm" | input$sh_plot == "rohs") {
        showNotification("please wait to run ROH", type = "message")
        print("please wait to run ROH")

        values$roh_mean <- roh(values$myphase, gwas = values$pieces,
                               w = input$wroh, fun = "mean",
                               snp_start = ival()[1], snp_end = ival()[2])
        values$roh_sum <- roh(values$myphase, gwas = values$pieces,
                              w = input$wroh, fun = "sum",
                              snp_start = ival()[1], snp_end = ival()[2])

        if (input$sh_plot == "rohm") p <- plot(values$roh_mean,
                                               title = "roh mean")
        if (input$sh_plot == "rohs") p <- plot(values$roh_sum,
                                               title = "roh_sum")
      }

      p
    }
  })




  # case_control choice for profile plot ------------------------------------


  output$res_choice <- renderUI({
    if (is.null(input$f_choice)) {
      return()
    }
    if (input$sh_plot != "profile") {
      return()
    }
    selectInput(
      "res_choice", NULL, choices = names(values$pieces$caco),
      selected = names(values$pieces$caco[1]), width = "300px")
  })



  # snp position  -----------------------------------------------------------

  output$pos <- renderUI({
    if (is.null(input$f_choice)) {
      return()
    }
    rng <- range(values$pieces[["snp.pos"]])
    if (sum(ival())) val <- ival() else val <- rng
    sliderInput("pos", "Position", min = rng[1], max = rng[2], value = val,
                step = 300, width = "900")
  })


  # show (write how many snp is selected) -----------------------------------

  output$snp_selected2 <- output$snp_selected <- renderText({
    len <- sum(values$pieces[["snp.pos"]] >= input$pos[1] &
                 values$pieces[["snp.pos"]] <= input$pos[2])
    if (input$sh_plot == "profile" && len >= 5000) {
      text <- "<b>Range of SNPs is longer than 5000 SNPs!</b>"
    } else {
      text <- paste("<b>Covered", len, "SNPs in this region</b>")
    }
    return(text)
  })




  output$mcnmr_choice <- renderUI({
    input$f_choice
    selectInput("mcnmr_choice", NULL, choices = names(values$caco),
                selected = names(values$caco)[1], width = "260px")
  })

  output$wmcn <- renderUI({
    sliderInput("wmcn", NULL, min = 1, max = 20, step = 1, value = 1,
                width = "140px")
  })


  output$nperm <- renderUI({
   if (!(substr(input$sh_plot, 1, 3) == "roh" &&
     input$roh_perm) || input$Panel != "nIBD") {
     return()
   }
   numericInput("nperm", NULL, value = 1000, min = 100, max = 100000,
                step = 1, width = "200px")
 })

  output$runtestr <- renderUI({
    input$wroh
    input$mcnmr_choice
    input$chr
    input$wmcn
    input$pos
    input$nperm
    if (!(substr(input$sh_plot, 1, 3) == "roh" &&
          input$roh_perm) || input$Panel != "nIBD") {
      return()
    }
    actionButton("runtestr", paste0("Test it!"))
  })

  runtestr <- eventReactive(input$runtestr, {
    myroh_mat <- roh(values$myphase, gwas = values$pieces, w = input$wroh,
                     fun = "sum", snp_start = ival()[1], snp_end = ival()[2],
                     roh_mat = T)
    mymcnemar <- mcnemar_test(values$roh_sum, reference = input$mcnmr_choice,
                              w = input$wmcn)
    mcnemar_test_permut(mymcnemar, roh_mat = myroh_mat, gwas = values$pieces,
                        nperm = input$nperm, reference = input$mcnmr_choice,
                        w = input$wmcn)
  })


  output$s_pos <- renderUI({
    rng <- range(values$pieces[["snp.pos"]])
    val <- ival()
    numericInput("s_pos", "Start Pos.", min = rng[1], max = rng[2],
                 value = val[1], width = "400")
  })

  output$e_pos <- renderUI({
    rng <- range(values$pieces[["snp.pos"]])
    val <- ival()
    numericInput("e_pos", "End Pos.", min = rng[1], max = rng[2],
                 value = val[2], width = "400")
  })

  output$wsize <- renderUI({
    max.win <- length(which(values$pieces[["snp.pos"]] >= ival()[1] &
                              values$pieces[["snp.pos"]] <= ival()[2])) - 1
    sliderInput("wsize", "Window Size", min = 2, max = min(100, max.win),
                step = 1, value = min(10, max.win), width = "200px")
  })

  observeEvent(input$pos, {
    ival(input$pos)
  })
  observeEvent(input$s_pos, {
    ival(c(input$s_pos, input$e_pos))
  })
  observeEvent(input$e_pos, {
    ival(c(input$s_pos, input$e_pos))
  })



  output$runit <- renderUI({
    actionButton("runit", paste0("RUN IBDW (Chr: ", input$chr, ") "))
  })



  observeEvent(input$runit, {
    showNotification("program is running ...", type = "message")
    print("program is running ...")
    values$hap_str <- list()
    values$hap_str <- haplotype_structure(values$myregion2,
                                          phase = values$myphase,
                                          w = input$wsize,
                                          snp_start = input$pos[1],
                                          snp_end = input$pos[2])
    values$happ_freq <- haplotype_frequency(values$hap_str)
    values$mygwid_win <- extract_window(values$myregion2, w = input$wsize,
                                        snp_start = input$pos[1],
                                        snp_end = input$pos[2])
    values$hap_str_genot <- list()
    values$hap_str_genot <- haplotype_structure(values$pieces,
                                                phase = values$myphase,
                                                w = input$wsize,
                                                snp_start = input$pos[1],
                                                snp_end = input$pos[2])
    values$happ_freq_genot <- haplotype_frequency(values$hap_str_genot)
  })

  output$IBD_out <- plotly::renderPlotly({
    if (is.null(values$hap_str)) {
      return()
    }
    if (input$sh_plot == "snp_win") {
      if (is.null(values$hap_str)) {
        return()
      }
      plot(values$mygwid_win, title = "gwid subset")
    } else if (input$sh_plot == "hap_win" & input$genome_ibd == FALSE) {
      plot(values$happ_freq, plot_type = "result_snps", title = "hap most freq")
    } else if (input$sh_plot == "hap_1" & input$genome_ibd == FALSE) {
      if (is.null(input$snp_win)) {
        return()
      }
      plot(values$happ_freq, plot_type = "haplotype_structure_frequency",
           nwin = input$snp_win, type = "version1")
    } else if (input$sh_plot == "hap_2" & input$genome_ibd == FALSE) {
      if (is.null(input$snp_win)) {
        return()
      }
      plot(values$happ_freq, plot_type = "haplotype_structure_frequency",
           nwin = input$snp_win, type = "version2")
    } else if (input$sh_plot == "hap_win" & input$genome_ibd == TRUE) {
      plot(values$happ_freq_genot, plot_type = "result_snps",
           title = "hap most freq")
    } else if (input$sh_plot == "hap_1" & input$genome_ibd == TRUE) {
      if (is.null(input$snp_win)) {
        return()
      }

      plot(values$happ_freq_genot, plot_type = "haplotype_structure_frequency",
           nwin = input$snp_win, type = "version1")
    } else if (input$sh_plot == "hap_2" & input$genome_ibd == TRUE) {
      if (is.null(input$snp_win)) {
        return()
      }

      plot(values$happ_freq_genot, plot_type = "haplotype_structure_frequency",
           nwin = input$snp_win, type = "version2")
    } else {
      return()
    }
  })


  output$snp_win <- renderUI({
    input$f_choice
    if (!input$sh_plot %in% c("hap_1", "hap_2")) {
      return()
    }
    idw <- length(which(values$pieces[["snp.pos"]] >= input$s_pos &
                          values$pieces[["snp.pos"]] <= input$e_pos))
    sliderInput("snp_win", "Select Window", min = 1,
                max = idw - input$wsize + 1, value = 1, step = 1, width = "900")
  })

  output$fisher_choice <- renderUI({
    if (!input$fisher) {
      return()
    }

    selectInput("fisher_choice", "Test", choices = names(values$pieces$caco),
                selected = names(values$pieces$caco)[1], width = "250px")
  })

  output$fisher_side <- renderUI({
    if (!input$fisher || substr(input$sh_plot, 1, 3) %in% c("hap", "sub")) {
      return()
    }
    selectInput("fisher_side", "Alt.",
                choices = c(">" = "greater", "!=" = "two.sided", "<" = "less"),
                selected = "greater", width = "150px")
  })
  output$nperm_gwid <- renderUI({
    if (!input$perm_gap) {
      return()
    }
    numericInput("nperm_gwid", NULL, value = 1000, min = 100, max = 100000,
                 step = 1, width = "200px")
  })

  output$runtest <- renderUI({
    if (input$Panel != "nIBD") {
      return()
    }

    if (input$fisher == TRUE) {
      actionButton("runtest", paste0("Test it! (",
                                     input$fisher_choice, "-",
                                     input$chr, "):", input$fisher_side, ":",
                                     input$pos[1] %% 1000, " - ",
                                     input$pos[2] %% 1000))
    } else if (input$fisher == FALSE & input$perm_gap == TRUE) {
      actionButton("runtest", paste0("Test it!"))
    }
  })

  output$runtestw <- renderUI({
    if (input$fisher == TRUE) {
      actionButton("runtestw", paste0("Win.Test (", input$fisher_choice, "-",
                                      input$chr, "):", input$fisher_side, ":",
                                      input$s_pos %% 1000, " - ",
                                      input$e_pos %% 1000))
    } else if (input$fisher == FALSE & input$perm_gap == TRUE) {
      actionButton("runtestw", paste0("Test it!"))
    }
  })

  output$download.b <- renderUI({
    fname <- paste("Chr", input$chr, ":", input$pos[1] %% 1000, "-",
                   input$pos[2] %% 1000)
    downloadButton("saveit", paste0("Downl. res. (", fname, ")"))
  })

  getReactiveData <- reactive({
    list(
      profile = values$myprofile,
      hap_str = values$hap_str,
      happ_freq = values$happ_freq,
      hap_str_genot = values$hap_str_genot,
      happ_freq_genot = values$happ_freq_genot,
      fish_test_gwid = values$fish_test_gwid,
      fish_test_gwas = values$fish_test_gwid,
      permutation_gwid = values$permutation_gwid,
      permutation_gwas = values$permutation_gwas,
      fisher_test_window_gwid = values$fisher_test_window_gwid
    )
  })

  output$saveit <- downloadHandler(
    filename = function() {
      "gwid_data_files.rda"
    },
    content = function(file) {
      data_to_save <- getReactiveData()
      saveRDS(data_to_save, file = file)
    }
  )

  runtest <- eventReactive(input$runtest, {
    showNotification("program is running ...", type = "message")
    print("program is running ...")

    if (input$fisher == TRUE) {
      if (input$sh_plot == "snp_ibd") {
        fishgwid <- fisher_test(values$myregion2, caco = values$caco,
                                snp_start = ival()[1], snp_end = ival()[2],
                                reference = input$fisher_choice,
                                alternative = input$fisher_side)
        values$fish_test_gwid <- fishgwid
        p <- plot(fishgwid, title = "gwid fisher test")
      } else if (input$sh_plot == "profile") {
        fishgwid <- fisher_test(values$myregion2, caco = values$caco,
                                snp_start = ival()[1], snp_end = ival()[2],
                                reference = input$fisher_choice,
                                alternative = input$fisher_side)
        values$fish_test_gwid <- fishgwid
        p <- plot(fishgwid, title = "gwid fisher test")
      } else if (input$sh_plot == "snp_cnt") {
        fishgwas <- fisher_test(values$pieces, reference = input$fisher_choice,
                                snp_start = ival()[1], snp_end = ival()[2],
                                alternative = input$fisher_side)
        values$fish_test_gwas <- fishgwas
        p <- plot(fishgwas, title = "gwas fisher test")
      }
    } else if (input$perm_gap == TRUE) {
      if (input$sh_plot == "snp_ibd" | input$sh_plot == "profile") {
        myperm_gap <- permutation_test(
          values$myregion2, values$pieces,
                                       snp_start = ival()[1],
                                       snp_end = ival()[2],
                                       nperm = input$nperm_gwid,
          reference = "cases", ibd_data = paste0(data_path(), "chr",
                                                 input$chr, ".ibd"))
        values$permutation_gwid <- myperm_gap
        p <- plot(myperm_gap, title = "permutation test gwid")
      }
      if (input$sh_plot == "snp_cnt") {
        myperm_gap <- permutation_test(
          values$pieces, snp_start = ival()[1], snp_end = ival()[2],
          nperm = input$nperm_gwid, reference = "cases")
        values$permutation_gwas <- myperm_gap
        p <- plot(myperm_gap, title = "permutation test gwas")
      }
    } else {
      return()
    }
    p
  })

  runtestw <- eventReactive(input$runtestw, {
    if (input$fisher == TRUE) {
      if (input$sh_plot == "snp_win") {
        fish_gwid_win <- fisher_test(
          values$mygwid_win, caco = values$caco,
          reference = input$fisher_choice, alternative = input$fisher_side)
        values$fisher_test_window_gwid <- fish_gwid_win
        p <- plot(fish_gwid_win, title = "gwid window fisher test")
      } else if (input$sh_plot == "hap_win" & input$genome_ibd == FALSE) {
        p <- plot(gtest(values$hap_str, reference = input$fisher_choice),
                  title = "gtest")
      } else if (input$sh_plot == "hap_1" & input$genome_ibd == FALSE) {
        p <- plot(gtest(values$hap_str, reference = input$fisher_choice),
                  title = "gtest")
      } else if (input$sh_plot == "hap_2" & input$genome_ibd == FALSE) {
        p <- plot(gtest(values$hap_str, reference = input$fisher_choice),
                  title = "gtest")
      } else if (input$sh_plot == "hap_win" & input$genome_ibd == TRUE) {
        p <- plot(gtest(values$hap_str_genot, reference = input$fisher_choice),
                  title = "gtest")
      } else if (input$sh_plot == "hap_1" & input$genome_ibd == TRUE) {
        p <- plot(gtest(values$hap_str_genot, reference = input$fisher_choice),
                  title = "gtest")
      } else if (input$sh_plot == "hap_2" & input$genome_ibd == TRUE) {
        p <- plot(gtest(values$hap_str_genot, reference = input$fisher_choice),
                  title = "gtest")
      } else {
        return()
      }
    }
    p
  })

  output$test_out <- plotly::renderPlotly({
    if (input$roh_test == TRUE) {
      p <- plot(mcnemar_test(values$roh_sum, reference = input$mcnmr_choice,
                             w = input$wmcn))
    } else if (input$roh_perm == TRUE) {
      p <- plot(runtestr())
    } else {
      runtest()
    }
  })

  output$test_outw <- renderPlotly({
    runtestw()
  })
}




#' laucnh a shiny app
#'
#' @param data_folder_address address of the folder that your data folders are.
#' for example if you have two sets of data such as data1 and data2 and they are
#' in mydata folder then your data_folder_address should be "./mydata"
#' @param ... other variables
#'
#' @return open a shiny app
#' @export
launch_app <- function(data_folder_address, ...) {
  shiny_env <- new.env()
  environment(shiny_ui) <- shiny_env
  environment(shiny_server) <- shiny_env
  assign("data_folder_address", data_folder_address, shiny_env)
  shiny::shinyApp(
    ui = shiny_ui,
    server = shiny_server
  )
}
