ui.modules_dash <- function(id) {
  ns <- NS(id)
  tagList(
    bs4Jumbotron(
      title = "Welcome to the Shiny TF-Target Finder (TFTF)!",
      lead = HTML("Transcription factors (TFs) are crucial in modulating gene expression and sculpting cellular and organismal phenotypes. The identification of TF-target gene interactions is pivotal for compre-hending molecular pathways and disease etiologies but has been hindered by the demanding nature of traditional experimental approaches. This paper introduces a novel web application and package, utilizing the R program, which predicts TF-target gene relationships and vice versa. Our application integrates the predictive power of various bioinformatic tools, leveraging their combined strengths to provide robust predictions. It merges databases for enhanced precision, incorporates gene expression correlation for accuracy, and employs pan-tissue correlation anal-ysis for context-specific insights. The application also enables the integration of user data with established resources to analyze TF-target gene networks.</br></br>Install PCAS packages: <b>remotes::install_github('WangJin93/TFTF')</b>, Then execute commands: <b>PCAS::TFTF_app() </b> to run the app locally.</br><b>Citation: </b>Wang J. TFTF: An R-Based Integrative Tool for Decoding Human Transcription Factor–Target Interactions. Biomolecules. 2024; 14(7):749. https://doi.org/10.3390/biom14070749"),
      status = "info",
      btnName = "Detailed introduction",
      href = "https://github.com/WangJin93/TFTF"
    ),
    fluidRow(
      column(6,
             bs4UserCard(
               title = bs4UserDescription(
                 title = h4("Application Developer"),
                 subtitle = "Jin wang",
                 image = "P1102513.jpg",
                 type = 1
               ),
               status = "info",
               width = 12,

               bs4ListGroup(
                 width = 12,
                 type = "action",
                 bs4ListGroupItem(
                   "Email: Jinwang93@suda.edu.cn"
                 ),
                 bs4ListGroupItem(
                   "Affiliation: Soochow University"
                 ),
                 bs4ListGroupItem(
                   "Personal website: https://www.jingege.wang",
                   href = "https://www.jingege.wang"
                 )
               )
             )
      ),
      column(6,
             box(title = h4("Citation"), solidHeader = T, status = "info", width = 12, collapsible = T,
                 h5("Wang J. TFTF: An R-Based Integrative Tool for Decoding Human Transcription Factor–Target Interactions. Biomolecules. 2024; 14(7):749. https://doi.org/10.3390/biom14070749")),
             box(title = h4("Package source code"), solidHeader = T, status = "info", width = 12, collapsible = T,
                 h5("https://github.com/WangJin93/TFTF"))
                 # h4("TF list obtained by taking the union of multiple tools"),
             # hr(),
             # shinycssloaders::withSpinner(DTOutput(outputId = ns("tf_info"))),
             #   downloadButton(ns("download_tf.csv"), "Download csv table",class = "mybutton")
             #
      )

    )
  )
}
server.modules_dash <- function(input, output, session) {
  ns <- session$ns

  output$tf_info <- renderDT({

    DT::datatable(
      tf_list,rownames = F,
      options = list(pageLength = 5)
    )

  })

  output$download_tf.csv <- downloadHandler(
    filename = function() {
      "TF_list.csv"
    },
    content = function(file) {
      write.csv(tf_list, file, row.names = FALSE)
    }
  )

}

