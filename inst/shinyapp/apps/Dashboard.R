ui.modules_dash <- function(id) {
  ns <- NS(id)
  tagList(
    bs4Jumbotron(
      title = "Welcome to the Shiny TF-Target Finder (TFTF)!",
      lead = "This is a Shiny application that bridging multiple predictive models for decoding transcription factor-target interactions in human. Our application synergizes the predictive capacities of multiple web tools by intersecting their results to enhance reliability!",
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
                 h5("Undering review")),
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

