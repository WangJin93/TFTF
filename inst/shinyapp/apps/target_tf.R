gene_list <- readRDS("rdata/pancan_gene_list.RDS")
ui.modules_target <- function(id) {
  ns <- NS(id)
  fluidPage(
    HTML("<h2 style='text-align: center;color: #666'>Transcription factor prediction</h2>

<hr style='border:3 double #987cb9' width='80%' color='#987cb9' size='10'>"),

    # Application title

    # Sidebar with a slider input for number of bins
    sidebarLayout(
      sidebarPanel(width = 3,

                   selectizeInput(
                     inputId = ns("Gene"),
                     label = "Input a gene symbol:",
                     choices = NULL,
                     width = "100%",
                     options = list(
                       create = TRUE
                     )
                   ),

                   selectizeInput(
                     inputId = ns("tables"),
                     label = "Select datasets:",
                     choices = c("hTFtarget",
                                 "KnockTF",
                                 "FIMO_JASPAR",
                                 "PWMEnrich_JASPAR",
                                 "ENCODE",
                                 "CHEA",
                                 "TRRUST",
                                 "GTRD",
                                 "ChIP_Atlas"),
                     multiple = T,
                     selected  = c("hTFtarget","KnockTF","FIMO_JASPAR"),
                     width = "100%"
                   ),
                   conditionalPanel(ns = ns,
                                    condition = "input.tables !== null && input.tables.indexOf('KnockTF') >= 0",
                                    sliderInput(ns("cut.log2FC"),"Log2FC Threshold for KnockTF:",min = 0,max = 3,value = 0.5,step = 0.5),
                                    checkboxInput(ns("down.only"),"Down-regulatted genes only?",F)
                   ),


                   hr(),
                   selectizeInput(
                     inputId = ns("correlation"),
                     label = "Correlaiton:",
                     choices = c("Corr_TCGA"= "TCGA",
                                 "Corr_GTEx" = "GTEx"),
                     multiple = T,
                     selected  = NULL,
                     width = "100%"
                   ),

                   uiOutput(ns("cor.setting")),
                   hr(),
                   shinyWidgets::actionBttn(
                     inputId = ns("search_bttn"),
                     label = "Go!",
                     style = "gradient",
                     icon = icon("search"),
                     color = "primary",
                     block = TRUE,
                     size = "sm"
                   ),
                   hr(),
                   selectizeInput(
                     inputId = ns("tables_venn"),
                     label = "Select datasets to get intersection:",
                     choices = c("hTFtarget",
                                 "KnockTF",
                                 "TRRUST",
                                 "ENCODE",
                                 "FIMO_JASPAR",
                                 "CHEA",
                                 "GTRD",
                                 "ChIP_Atlas"),
                     multiple = T,
                     selected  = c(""),
                     width = "100%"
                   ),

      ),
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(
          tabPanel('All results',
                   br(),
                   shinyjs::useShinyjs(),  # Set up shinyjs
                   HTML("<hr>"),
                   shinycssloaders::withSpinner(DTOutput(outputId = ns("results"))),
                   shinyjs::hidden(
                     shinyWidgets::downloadBttn(ns("download.csv"), "Download csv table")
                   )),
          tabPanel('Venn diagram',
                   br(),
                   shinyjs::useShinyjs(),  # Set up shinyjs
                   HTML("<hr>"),
                   fluidRow(
                     column(4,
                          sliderInput(ns('linewidth'), 'Select line width', min = 0, max = 5, value = 1, step = 0.1),
                          radioButtons(ns("linetype"), "Select line type:", choices = c("solid" = 1,"dashed" = 2, "dotted" = 3,

                                                                                    "dot-dashed" = 4, "longdashed" = 5, "blank" = 0), inline = TRUE),
                          selectizeInput(
                            inputId = ns("color_panel"),
                            label = "Select color_panel for plotting:",
                            choices = c("Set1",
                                        "Set2",
                                        "Set3",
                                        "Accent",
                                        "Dark2",
                                        "Paired",
                                        "Pastel1",
                                        "Pastel2"),
                            multiple = F,
                            selected  = c("Set1"),
                            width = "100%"
                          ),
                          sliderInput(ns('fontsize'), 'Select font size', min = 1, max = 2, value = 1, step = 0.1),
                   ),
                   column(8,align="center",
                          shinycssloaders::withSpinner(plotOutput(ns("venn_diagram"), height = "600px", width = "600px")),
                          shinyjs::hidden(
                            shinyWidgets::downloadBttn(ns("download"), "Download Figure")
                          )

                   )#column
                   )
          ),
          tabPanel('Individual dataset',
                   br(),
                   HTML("<hr>"),
                   selectizeInput(
                     inputId = ns("individual"),
                     label = "Select dataset:",
                     choices = c("hTFtarget",
                                 "KnockTF",
                                 "TRRUST",
                                 "ENCODE",
                                 "FIMO_JASPAR",
                                 "CHEA",
                                 "GTRD",
                                 "ChIP_Atlas"),
                     multiple = F,
                     selected  = "hTFtarget",
                     width = "100%",
                   ),

                   HTML("<hr>"),
                   shinycssloaders::withSpinner(DTOutput(outputId = ns("individual_data"))),
                   hr(),
                   shinyjs::hidden(
                     shinyWidgets::downloadBttn(ns("download.individual"), "Download individual data")
                   ),
                   tags$head(tags$style(".mybutton{background-color:aliceblue;} .mybutton2{background-color:antiquewhite;} .skin-black .sidebar .mybutton{color: green;}") )

          )
        )
      )
    )
  )
  }
server.modules_target <- function(input, output, session) {
  ns <- session$ns
  w <- waiter::Waiter$new(id = ns("gene_pancan_dist"), html = waiter::spin_hexdots(), color = "black")
  observe({
    updateSelectizeInput(
      session,
      "Gene",
      choices = gene_list,
      selected = "GAPDH",
      server = TRUE
    )
  })
  observe({
    if (length(input$correlation) == 2){
      output$cor.setting <- renderUI({
        tagList(

          selectizeInput(
            inputId = ns("TCGA_type"),
            label = "TCGA cancer type:",
            choices = tissue$TCGA,
            multiple = F,
            selected  = "LUAD",
            width = "100%"
          ),

          selectizeInput(
            inputId = ns("GTEx_type"),
            label = "GTEx tissue type:",
            choices = tissue$GTEx,
            multiple = F,
            selected  = "Lung",
            width = "100%"

          ),
          numericInput(ns("cor.threshold"),
                       "coefficient threshold",
                       value = 0.3,
                       min = 0.3,
                       max = 1,
                       step = 0.1)
        )
      }
      )}else if (length(input$correlation) == 1){
        if (input$correlation == "TCGA"){
          output$cor.setting <- renderUI({
            tagList(

              selectizeInput(
                inputId = ns("TCGA_type"),
                label = "TCGA cancer type:",
                choices = tissue$TCGA,
                multiple = T,
                selected  = "LUAD",
                width = "100%"
              ),
              numericInput(ns("cor.threshold"),
                           "coefficient threshold",
                           value = 0.3,
                           min = 0.3,
                           max = 1,
                           step = 0.1)
            )
          })
        }
        if(input$correlation == "GTEx"){
          output$cor.setting <- renderUI({
            tagList(
              selectizeInput(
                inputId = ns("GTEx_type"),
                label = "GTEx tissue type:",
                choices = tissue$GTEx,
                multiple = T,
                selected  = "Lung",
                width = "100%"
              ),
              numericInput(ns("cor.threshold"),
                           "coefficient threshold",
                           value = 0.3,
                           min = 0.3,
                           max = 1,
                           step = 0.1)
            )
          })
        }}else{
          output$cor.setting <- renderUI({
            tagList(
              NULL)
          })

        }
  })

  TF_results <- eventReactive(input$search_bttn,{
    TF_results <- predict_TF(datasets = input$tables,
                                   target = input$Gene,
                                   cut.log2FC = input$cut.log2FC,
                                   down.only = input$down.only,
                                   cor_DB = input$correlation,
                                   TCGA_tissue = input$TCGA_type,
                                   GTEx_tissue = input$GTEx_type,
                                   cor_cutoff = input$cor.threshold,
                             app=T)

    return(TF_results)
  })
  observe({
    updateSelectInput(session, "individua",
                      choices =  names(TF_results()$results),
                      selected =  names(TF_results()$results)[1]
    )
    updateSelectInput(session, "tables_venn",
                      choices =  names(TF_results()$results),
                      selected = names(TF_results()$results)
    )

  })

  observeEvent(input$search_bttn, {
    shinyjs::show(id = "download")
    shinyjs::show(id = "download.csv")
    shinyjs::show(id = "download.individual")

  })
  output$results <- renderDT({
    tf_results <- TF_results()
    tf_results <- intersections(tf_results,datasets_index = which(names(tf_results) %in% input$tables_venn))
    for (i in names(tf_results)) {
      length(tf_results[[i]]) <- lapply(tf_results, length) %>% unlist() %>% max()
    }

    tf_results_d <- as.data.frame(tf_results)
    DT::datatable(
      tf_results_d
    )

  })

  output$venn_diagram <- renderPlot({
    tf_results <- TF_results()
    inter_data <- intersections(tf_results,datasets_index = which(names(tf_results) %in% input$tables_venn))
    dd <- which(names(inter_data)=="intersection")
    if (length(dd)!=0){
      inter_data <- inter_data[-which(names(inter_data)=="intersection")]
    }
    if (length(inter_data) < 6){
      VD <- venn.diagram(inter_data, filename = NULL, fill = RColorBrewer::brewer.pal(length(inter_data), input$color_panel), cex = input$fontsize,
                         margin = 0.2, cat.cex = input$fontsize, lwd = input$linewidth,
                         lty = rep(as.numeric(input$linetype), length(inter_data)))
      grid.draw(VD)

    }else{
      flowerplot(inter_data,
                 ellipse_col_pal = input$color_panel,
                 circle_col = "white",
                 label_text_cex = input$fontsize)


    }


  })
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$Gene,"_TFs.pdf")
    },
    content = function(file) {

        pdf(file)
      tf_results <- TF_results()
      inter_data <- intersections(tf_results,datasets_index = which(names(tf_results) %in% input$tables_venn))
      dd <- which(names(inter_data)=="intersection")
      if (length(dd)!=0){
        inter_data <- inter_data[-which(names(inter_data)=="intersection")]
      }
      if (length(inter_data) < 6){
        VD <- venn.diagram(inter_data, filename = NULL, fill = RColorBrewer::brewer.pal(length(inter_data), input$color_panel), cex = input$fontsize,
                           margin = 0.2, cat.cex = input$fontsize, lwd = input$linewidth,
                           lty = rep(as.numeric(input$linetype), length(inter_data)))
        grid.draw(VD)

      }else{
        flowerplot(inter_data,
                   ellipse_col_pal = input$color_panel,
                   circle_col = "white",
                   label_text_cex = input$fontsize)


      }



      dev.off()

    })

  output$download.csv <- downloadHandler(
    filename = function() {
      paste0(input$Gene, "_TF_prediction.csv")
    },
    content = function(file) {
      tf_results <- TF_results()
      tf_results <- intersections(tf_results,datasets_index = which(names(tf_results) %in% input$tables_venn))
      for (i in names(tf_results)) {
        length(tf_results[[i]]) <- lapply(tf_results, length) %>% unlist() %>% max()
      }

      tf_results_d <- as.data.frame(tf_results)
      write.csv(tf_results_d, file, row.names = FALSE)
    }
  )

  output$individual_data <- renderDT({
    tf_results <- TF_results()
    DT::datatable(
      tf_results[[input$individual]]
    )
  })

  output$download.individual <- downloadHandler(
    filename = function() {
      paste0(input$Gene, "_TF_",input$individual,"_prediction.csv")
    },
    content = function(file) {
      tf_results <- TF_results()
      write.csv(tf_results[[input$individual]], file, row.names = FALSE)
    }
  )
}
