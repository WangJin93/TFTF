ui.modules_tf_list <- function(id) {
  ns <- NS(id)
  mainPanel(width = 10,
    h4("We only include transcription factors that are present in at least 2 datasets out of nine, resulting in a total of 1575 transcription factors."),
    hr(),
    shinycssloaders::withSpinner(DTOutput(outputId = ns("tf_info"))),
    br(),
  shinyWidgets::downloadBttn(ns("download_tf.csv"), "Download csv table")
)
  }
server.modules_tf_list <- function(input, output, session) {
  ns <- session$ns

  output$tf_info <- renderDT({

   datatable(
      tf_list,rownames = F
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

