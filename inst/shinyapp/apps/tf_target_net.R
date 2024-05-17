
filetype_map <- c( 'csv' = ',', 'txt'="\t")
hTFtarget <- tf_list %>% dplyr::filter(hTFtarget == "T") %>% .[,1]
KnockTF <- tf_list %>% dplyr::filter(.,KnockTF == "T") %>% .[,1]
FIMO_JASPAR <- tf_list %>% dplyr::filter(.,FIMO_JASPAR == "T") %>% .[,1]
PWMEnrich_JASPAR  <- tf_list %>% dplyr::filter(.,PWMEnrich_JASPAR  == "T") %>% .[,1]
CHEA <- tf_list %>% dplyr::filter(.,CHEA == "T") %>% .[,1]
ENCODE <- tf_list %>% dplyr::filter(.,ENCODE == "T") %>% .[,1]
TRRUST <- tf_list %>% dplyr::filter(.,TRRUST == "T") %>% .[,1]
GTRD <- tf_list %>% dplyr::filter(.,GTRD == "T") %>% .[,1]
ChIP_Atlas <- tf_list %>% dplyr::filter(.,ChIP_Atlas == "T") %>% .[,1]

ui.modules_net <- function(id) {
  ns <- NS(id)
  fluidPage(
    HTML("<h2 style='text-align: center;color: #666'>TF-targets network</h2>

<hr style='border:3 double #987cb9' width='80%' color='#987cb9' size='10'>"),

    # Application title

    # Sidebar with a slider input for number of bins
    sidebarLayout(
      sidebarPanel(width = 3,
                   fileInput(ns('datafile'), "Choose text/csv File",
                             accept=c('text/csv',
                                      'text/comma-separated-values,text/plain',
                                      '.csv'),multiple = FALSE),
                   p("Please make sure your data contains following columns: 'logFC','P.Value','Symbol'"),
                   HTML("<a href='DEG_results.csv'> <i class='fa fa-download'> </i> Download example data</a> <hr>  "),
                   numericInput(ns("logFC_cutoff"), "|log2FC| threshold:", value = 1,min = 0),
                   selectInput(
                     ns("p_cutoff"),
                     "P value threshold:",
                     c("0.05","0.01","0.001","0.0001"),
                     selected = "0.05",
                     multiple = FALSE,
                     selectize = TRUE,
                     width = NULL,
                     size = NULL
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
                     selected  = c("hTFtarget","KnockTF"),
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
                   selectizeInput(
                     ns("tf_interest"),
                     "TF to analysis:",
                     choices = NULL,
                     multiple = TRUE,
                   ),



                   hr(),
                   shinyWidgets::actionBttn(
                     inputId = ns("search_bttn"),
                     label = "Go!",
                     style = "gradient",
                     icon = icon("search"),
                     color = "primary",
                     block = TRUE,
                     size = "sm"
                   )

      ),
      # Show a plot of the generated distribution
      mainPanel(
        bs4Dash::tabsetPanel(	id = ns("tablist"),
          tabPanel('TF results',
                   br(),
                   shinyjs::useShinyjs(),  # Set up shinyjs
                   shinycssloaders::withSpinner(DTOutput(outputId = ns("tf_results"))),
                   shinyWidgets::downloadBttn(ns("download_tf_results.csv"), "Download csv table")
                   ),
          tabPanel('TF-target network',	value = "net",
                   br(),
                   shinyjs::useShinyjs(),  # Set up shinyjs
                   fluidRow(
                     column(4,
                            selectInput(inputId = ns("layout"),
                                        label = "Layout style:",
                                        c("fr",
                                          "dh",
                                          "gem",
                                          "star",
                                          "lgl",
                                          "randomly",
                                          "mds",
                                          "grid",
                                          "kk",
                                          "graphopt",
                                          "drl",
                                          "circle"),
                                        selected = "fr"),


                            sliderInput(ns("text_size"), "Text size:",
                                        min = 0,
                                        max = 8,
                                        value = 3)
                            ),
                   column(8,

                          shinycssloaders::withSpinner(plotOutput(ns("network"), height = "auto", width = "600px")),
                            shinyWidgets::downloadBttn(ns("download"), "Download Figure"),
                          tags$head(tags$style(".mybutton{background-color:aliceblue;} .mybutton2{background-color:antiquewhite;} .skin-black .sidebar .mybutton{color: green;}") )

                   )#column
                   )
          ),
          tabPanel('Plotting data',
                   br(),

                   shinycssloaders::withSpinner(DTOutput(outputId = ns("plot.data"))),
                   hr(),
                     shinyWidgets::downloadBttn(ns("download.data"), "Download individual data"),
                   tags$head(tags$style(".mybutton{background-color:aliceblue;} .mybutton2{background-color:antiquewhite;} .skin-black .sidebar .mybutton{color: green;}") )

          )
        )
      )
    )
  )
}
server.modules_net <- function(input, output, session) {
  ns <- session$ns
  w <- waiter::Waiter$new(id = ns("TF"), html = waiter::spin_hexdots(), color = "black")
  volcano_data <- reactive({
    inFile <- input$datafile
    if (is.null(inFile)) {
      return(NULL)
    }else{
      file_extension =  unlist(strsplit(inFile$datapath, '[.]'))[length(unlist(strsplit(inFile$datapath, '[.]')))]
      if(file_extension %in% names(filetype_map)){
        filetype <- filetype_map[file_extension]
        names(filetype) <- NULL

      } else {
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          paste("wrong file format", file_extension)
        ))
        return(NULL)
      }
      expr <-read.delim(inFile$datapath, header = TRUE, sep = filetype)

    }
    if ( all(c("logFC","P.Value","Symbol") %in% colnames(expr)) ){
      degs <- expr %>% dplyr::filter(abs(logFC)> input$logFC_cutoff & P.Value < as.numeric( input$p_cutoff) )

    }else{
      showModal(modalDialog(
        title = "Message",   easyClose = TRUE,
        "Please confirm the column names is corrected! 'Symbol' for gene column, 'logFC' for log2FC column and 'P.Value' for P value "
      ))
      return(NULL)

    }

    return(degs)
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
                multiple = F,
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
          })
        }}else{
          output$cor.setting <- renderUI({
            tagList(
              NULL)
          })

        }
  })

  gene_tf <- eventReactive({input$tables
    input$datafile
    input$p_cutoff
    input$logFC_cutoff},{
    datasets <- input$tables
    if (is.null(volcano_data())){
      return(NULL)
    }else{
      tf_DEG <- volcano_data() %>% dplyr::filter(Symbol %in% tf_list$TF)
      intersect( Reduce(intersect,lapply(datasets,function(x)  get(x))),tf_DEG$Symbol)
    }

  })
  observe({

    updateSelectizeInput(
      session,
      "tf_interest",
      choices = gene_tf(),
      selected = gene_tf()[1:2]
      )
  })


  output$tf_results <- renderDT({
    if (is.null(volcano_data())){NULL}else{
      DT::datatable(
        volcano_data() %>% dplyr::filter(Symbol %in% tf_list$TF)
      )
    }

  })
  observeEvent(input$search_bttn, {
    updateTabsetPanel(session, inputId = "tablist", selected = "net")
  })

  plot_data <- eventReactive(input$search_bttn, {
    if (is.null(volcano_data())) {
      return(NULL)
    }else{
      datasets <- input$tables
      all_results <- TF_Target_batch(datasets,tfs = input$tf_interest,
                                     cut.log2FC = input$cut.log2FC,
                                     down.only = input$down.only,
                                     cor_DB = input$correlation,
                                     TCGA_tissue = input$TCGA_type,
                                     GTEx_tissue = input$GTEx_type,
                                     cor_cutoff = input$cor.threshold,
                                     app = T)
    }
    all_results <- all_results %>% dplyr::filter(Target %in% volcano_data()$Symbol)
    return(all_results)
  })

  output$network <- renderPlot(height = 600,{

    all_results <- na.omit(plot_data())
    colnames(all_results) <- c("from", "to")
    graph_gt <- as_tbl_graph(all_results)
    nodes <- graph_gt %>% as.data.frame() %>% merge(.,volcano_data()[c("Symbol","logFC")] %>% unique(),by.x="name",by.y = "Symbol",all.x=T)
    graph_gt <- graph_gt  %>%
      mutate(logFC = nodes$logFC)


    ggraph(graph_gt, layout = input$layout) +
      geom_edge_fan(color="grey",show.legend=FALSE) +
      geom_node_point(aes(size = abs(logFC),fill = logFC),colour="black",shape=21)+
      geom_node_text(aes(label=name),size=input$text_size,repel=TRUE)+
      theme_graph() +guides(size="none")+
      scale_fill_distiller(palette = "Spectral")+
      theme(text = element_text(family = "Times"))

  })
  output$download <- downloadHandler(
    filename = function() {
      paste0("TF-target network.pdf")
    },

    content = function(file) {
      all_results <- na.omit(plot_data())
      colnames(all_results) <- c("from", "to")
      graph_gt <- as_tbl_graph(all_results)
      nodes <- graph_gt %>% as.data.frame() %>% merge(.,volcano_data()[c("Symbol","logFC")] %>% unique(),by.x="name",by.y = "Symbol", all.x = T)
      graph_gt <- graph_gt  %>%
        mutate(logFC = nodes$logFC)
      p <- ggraph(graph_gt, layout = input$layout) +
        geom_edge_fan(color="grey",show.legend=FALSE) +
        geom_node_point(aes(size = abs(logFC),fill = logFC),colour="black",shape=21)+
        geom_node_text(aes(label=name),size=input$text_size,repel=TRUE,max.overlaps=10)+
        theme_graph() +guides(size="none")+
        scale_fill_distiller(palette = "Spectral")+
        theme(text = element_text(family = "Times"))



      pdf(file, width = 8 ,height = 8)
      print(p)
      dev.off()
    })

  output$download_tf_results.csv <- downloadHandler(
    filename = function() {
      "TF_DEGs.csv"
    },
    content = function(file) {
         write.csv(volcano_data() %>% dplyr::filter(Symbol %in% tf_list$TF),
                   file, row.names = FALSE)
    }
  )

  output$plot.data <- renderDT({
    DT::datatable(
      plot_data() %>% merge(.,volcano_data()[c("Symbol","logFC")] %>% unique(),by.x="Target",by.y = "Symbol", all.x = T)
    )
  })

  output$download.data <- downloadHandler(
    filename = function() {
      "TF-target network.csv"
    },
    content = function(file) {
      results <-       plot_data() %>% merge(.,volcano_data()[c("Symbol","logFC")] %>% unique(),by.x="Target",by.y = "Symbol", all.x = T)

      write.csv(results, file, row.names = FALSE)
    }
  )
}
