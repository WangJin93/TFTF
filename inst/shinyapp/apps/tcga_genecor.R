

ui.modules_tcga <- function(id) {
  ns <- NS(id)
  fluidPage(
    HTML("<h2 style='text-align: center;color: #666'>Pan-tissue correlation analysis</h2>

<hr style='border:3 double #987cb9' width='80%' color='#987cb9' size='10'>"),
    shinyjs::useShinyjs(),
    use_waiter(),
    fluidRow(
      column(
        3,
        wellPanel(
          selectizeInput(
            inputId = ns("Pancan_search1"),
            label = "Input a TF name",
            choices = NULL,
            width = "100%",
            options = list(
              create = F
            )
          ),

          selectizeInput(
            inputId = ns("Pancan_search2"),
            label = "Input a gene symbol",
            choices = NULL,
            width = "100%",
            options = list(
              create = F,
              maxOptions = 5,
              placeholder = "Enter a gene symbol, e.g. TP53",
              plugins = list("restore_on_backspace")
            )
          ),
          selectInput(
            inputId = ns("data_source"),
            label = "Select database",
            choices = c("TCGA","GTEx","CCLE"),
            selected = "TCGA", multiple = F
          ),
          conditionalPanel(ns=ns,
                           condition = "input.data_source == 'TCGA'",
                           selectInput(
                             inputId = ns("data_type"),
                             label = "Only use normal or tumor",
                             choices = c("normal","tumor"),
                             selected = "tumor", multiple = TRUE
                           )
          ),
          selectInput(
            inputId = ns("cor_method"),
            label = "Select Correlation method",
            choices = c("spearman", "pearson"),
            selected = "pearson"
          ),
          shinyWidgets::sliderTextInput(
            inputId = ns("xinter"),
            label = "Coefficient cutpoint:",
            choices = seq(
              from = 0,
              to = 1,
              by = 0.05
            ),
            selected = 0.3,
            grid = TRUE
          ),
          selectInput(
            inputId = ns("yinter"),
            label = "P value cutpoint:",
            choices = c("0.05", "0.005","0.001"),
            selected = "0.05"
          ),
          br(),
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

        )),
      column(
        3,wellPanel(

          shinyWidgets::sliderTextInput(
            inputId = ns("alpha"),
            label = "Choose a transparent value",
            choices = seq(
              from = 0,
              to = 1,
              by = 0.1
            ),
            selected = "0.5",
            grid = TRUE
          ),
          shinyWidgets::sliderTextInput(
            inputId = ns("size"),
            label = "Dot size value",
            choices = seq(
              from = 0,
              to = 6,
              by = 0.5
            ),
            selected = 3,
            grid = TRUE
          ),

          textInput(ns('venncolors'), "Comma seperated list of colors", value = 'green,black,red'),
          a(href = "http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf",target = "_blank", "R-colors"),
          tags$hr(style = "border:none; border-top:2px solid #5E81AC;"),
        ),
        br(),
        br(),
        wellPanel(

          shinyWidgets::downloadBttn(
            outputId = ns("download"),
            style = "gradient",
            color = "default",
            block = TRUE,
            size = "sm"
          )
        )
      ),
      column(
        6,
        plotOutput(ns("gene_cor"), height = "auto"),
        hr(),
        h5("NOTEs:"),
        p("1. The data query may take some time based on your network. Wait until a plot shows"),
        p("2. ", tags$a(href = "https://pancanatlas.xenahubs.net/",target = "_blank", "Genomic profile data source")),
        p("3. Click the following data to view single cancer type "),
        tags$br(),
        DT::DTOutput(outputId = ns("tbl")),
          wellPanel(
            id = ns("save_csv"),
            shinyWidgets::downloadBttn(ns("downloadTable"), "Save as csv")
          )
      )
    )
  )
}
pancan_gene_list <- readRDS("rdata/pancan_gene_list.RDS")

server.modules_tcga <- function(input, output, session) {
  ns <- session$ns
  observe({
    updateSelectizeInput(
      session,
      "Pancan_search2",
      choices = pancan_gene_list,
      selected = "GAPDH",
      server = TRUE
    )
  })


  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("gene_cor"), html = waiter::spin_hexdots(), color = "white")

  # reformat input color
  plot.color <- reactive({
    venncolors <- input$venncolors
    color.tmp <- unlist(strsplit(venncolors, ','))
  })
  observe({
    updateSelectizeInput(
      session,
      "Pancan_search1",
      choices = tf_list$TF,
      selected = "STAT3",
      server = TRUE
    )
  })

  plot_data <- eventReactive(input$search_bttn, {
    if (nchar(input$Pancan_search1) >= 1 & nchar(input$Pancan_search2) >= 1) {
      output <- pantissue_cor_analysis(
        Gene1 = input$Pancan_search1,
        Gene2 = input$Pancan_search2,
        type = input$data_type,
        cor_method = input$cor_method,
        data_source = input$data_source
      )
    }
    if (is.null(output)){
      return(NULL)
    }
    return(output)
  })
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste0(input$Pancan_search1, "_", input$profile1, "_", input$Pancan_search2, "_", input$profile2, "_cor_",input$data_source,".csv")
    },
    content = function(file) {
      write.csv(plot_data()$cor_result,file,  row.names = FALSE)
    }
  )


  output$gene_cor <- renderPlot(  height = 600,{
    w$show() # Waiter add-ins
    output <- plot_data()$cor_result
    output$sigmark <- ifelse(output$p>0.05,"No",ifelse(output$r> input$xinter,"Positive",ifelse(output$r< -input$xinter,"Negative","No")))

    output$label=ifelse(output$sigmark != "No",as.vector(output$tissue),"")

    if (!is.null( output)){
      p <- ggplot(data = output,aes(x = r,y = logP)) +
        geom_point(alpha=input$alpha, size= input$size,aes(color=sigmark)) +
        scale_color_manual(values= plot.color())+      ##########color
        geom_vline(xintercept=c(input$xinter,-input$xinter),lty=4,col="black",lwd=0.8) +
        geom_hline(yintercept = -log10(as.numeric(input$yinter)) ,lty=4,col="black",lwd=0.8) +
        ylab("-log10(P value)")+xlab("Correlation coefficient")+
        theme_bw()

      p<- p + ggrepel::geom_label_repel(aes(label = label),data = output,    color="black"  )+
        theme(axis.title.x =element_text(size=14,    color="black" ), axis.title.y=element_text(size=14,    color="black" ),
              axis.text.x = element_text(size=14,    color="black" ), axis.text.y=element_text(size=14,    color="black"),
              legend.text=element_text(size=14),legend.title=element_text(size=14))


      p+ ggplot2::theme(
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15)
      )
    }

  })

  # download module
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$Pancan_search1, "_", input$profile1, "_", input$Pancan_search2, "_", input$profile2, "_cor_",input$data_source,".pdf")
    },
    content = function(file) {
      output <- plot_data()$cor_result
      output$sigmark <- ifelse(output$p>0.05,"No",ifelse(output$r> input$xinter,"Positive",ifelse(output$r< -input$xinter,"Negative","No")))

      output$label=ifelse(output$sigmark != "No",as.vector(output$tissue),"")
      p <- ggplot(data = output,aes(x = r,y = logP)) +
        geom_point(alpha=input$alpha, size= input$size,aes(color=sigmark)) +
        scale_color_manual(values= plot.color())+      ##########color
        geom_vline(xintercept=c(input$xinter,-input$xinter),lty=4,col="black",lwd=0.8) +
        geom_hline(yintercept = -log10(as.numeric(input$yinter)) ,lty=4,col="black",lwd=0.8) +
        ylab("-log10(P value)")+xlab("Correlation coefficient")+
        theme_bw()
      output$label=ifelse(output$sigmark != "No",as.vector(output$tissue),"")
      p<- p + ggrepel::geom_label_repel(aes(label = label),data = output,    color="black"  )+
        theme(axis.title.x =element_text(size=14,    color="black" ), axis.title.y=element_text(size=14,    color="black" ),
              axis.text.x = element_text(size=14,    color="black" ), axis.text.y=element_text(size=14,    color="black") )


      p <- p+ ggplot2::theme(
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15)
      )
      pdf(file,width = 7.5,height = 6)
      print(p)
      dev.off()
    }
  )
    # return data


    output$tbl <- DT::renderDataTable(
      plot_data()$cor_result, server = TRUE,selection = 'single'
    )

    # Keep selected database

    corr_func <- eventReactive(input$tbl_rows_selected, {
      s <- input$tbl_rows_selected
      df <- plot_data()
      if (length(s)) {
        selected <- df$cor_result[s,]
        #data<- df[which(df$tissue == type),]
        type <- selected$tissue
        print(type)
        df<-df$cor_data
        if (input$data_source == "TCGA"){
          df %>%
            dplyr::filter(.data$type1 == input$data_source) %>%
            dplyr::filter(.data$type2 %in% input$data_type) -> df
        }
        if (input$data_source == "GTEx"){
          df %>%
            dplyr::filter(.data$type1 == input$data_source) -> df
        }

        df <- df[which(df$tissue == type),]
        df <- df[c("sample","tissue",input$Pancan_search1,input$Pancan_search2)]
      }
      return(df)
    })



    observeEvent(eventExpr = input$tbl_rows_selected,{
      message(input$Pancan_search, " is queried by user from home search box.")
      s <- input$tbl_rows_selected

      if (length(s)) {
        showModal(
          modalDialog(
            title = paste("Gene-Gene correlation scatter plot"),
            size = "l",
            fluidPage(
              div(style="display:flex",
                  numericInput(inputId = ns("height_scatter"), label = "Height", value = 300,max = 600,width = 100),
                  numericInput(inputId = ns("width_scatter"), label = "Width", value = 500,width = 100),
                  textInput(inputId = ns('x_lab'), "X-axis suffix", value = ' relative expression',width = 250),
                  textInput(inputId = ns('y_lab'), "Y-axis suffix", value = ' relative expression',width = 250)

              ),
              column(
                12,
                plotOutput(ns("scatter_plot"),width = "100%",height = "auto",),
                DT::DTOutput(outputId = ns("scatter_corr")),
                wellPanel(
                  id = ns("save_scatter_csv"),
                  downloadButton(ns("downloadTable2"), "Save as csv")
                )

              ),
              column(
                12,
                h4("NOTEs:"),
                h5("1. The data query may take some time based on your network. Wait until a plot shows"),
                h5("2. The unit of gene expression is log2(tpm+0.001)"),
              )
            )
          )
        )
        width_scatter <- reactive ({ input$width_scatter })
        height_scatter <- reactive ({ input$height_scatter })
        output$scatter_plot <- renderPlot(width = width_scatter,
                                          height = height_scatter,{
                                            df<-corr_func()
                                            ggcorplot(df,input$Pancan_search1,input$Pancan_search2,
                                                      method=input$cor_method,x_lab=input$x_lab,y_lab=input$y_lab)
                                          })
        output$scatter_corr <- DT::renderDataTable(
          corr_func(), server = TRUE,selection = 'single'
        )
        ## downloadTable
        output$downloadTable2 <- downloadHandler(
          filename = function() {
            paste0(input$Pancan_search1, "_", input$Pancan_search2, "_", corr_func()[1,2], "_cor.csv")
          },
          content = function(file) {
            write.csv(corr_func(), file, row.names = FALSE)
          }
        )
      }
    })
}
