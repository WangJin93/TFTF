#' feedback UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_feedback_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             HTML("<h2 style='text-align: center;color: #666'>Feedback</h2>

<hr style='border:3 double #987cb9' width='80%' color='#987cb9' size='10'>"),
        bs4Card(
          title = "Details",
          status = "primary",
          solidHeader = FALSE,
          collapsible = FALSE,
          collapsed = FALSE,
          closable = FALSE,
          label = NULL,
          width = 12,
          tagList(
            fluidRow(
              column(6,
                textInput(
                  ns("issue_title"),
                  "Title",
                  value = "",
                  placeholder = "enter short label",
                  width = "100%"
                )
              ),
              column(6,
                     shinyWidgets::pickerInput(
                  ns("issue_labels"),
                  label = "Choose one category",
                  choices = c("bug", "enhancement"),
                  multiple = FALSE
                )
              )
            ),
            fluidRow(
              column(12,
                textAreaInput(
                  ns("issue_description"),
                  label = "Enter description",
                  value = "",
                  width = "300%",
                  cols = 80,
                  rows = 10,
                  placeholder = "Any suggestion or bugs feedback",
                  resize = "vertical"
                )
              )
            ),
            fluidRow(
              column(2,
                shinyWidgets::actionBttn(
                  ns("submit_issue"),
                  "Submit!",
                  icon = icon("save"),
                  style = "jelly",
                  color = "success",
                  size = "md"
                )
              )
            )
          )
        )
      )
    )
  )
}
    
#' feedback Server Function
#'
#' @noRd 
mod_feedback_server <- function(input, output, session){
  ns <- session$ns
  
  observeEvent(input$submit_issue, {
    # perform checks for mandatory inputs
    if (!shiny::isTruthy(input$issue_title)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Oops!",
        text = "Please enter a title for your issue.",
        type = "error"
      )
      return(NULL)
    }
    
    if (!shiny::isTruthy(input$issue_labels)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Oops!",
        text = "Please choose at least one category for your issue.",
        type = "error"
      )
      return(NULL)
    }
    
    if (!shiny::isTruthy(input$issue_description)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Oops!",
        text = "Please enter a description for your issue.",
        type = "error"
      )
      return(NULL)
    }
    
    # submit issue
    # Security note: never embed email credentials in source code. The SMTP
    # account/password are read from environment variables (or options) and
    # must be configured on the deployment machine, e.g.:
    #   Sys.setenv(TFTF_MAIL_FROM = "you@qq.com",
    #              TFTF_MAIL_USER = "you@qq.com",
    #              TFTF_MAIL_PASS = "<smtp authorization code>",
    #              TFTF_MAIL_TO   = "admin@example.com")
    mail_pass <- Sys.getenv("TFTF_MAIL_PASS", unset = getOption("TFTF.mail_pass", ""))
    mail_user <- Sys.getenv("TFTF_MAIL_USER", unset = getOption("TFTF.mail_user", ""))

    if (!nzchar(mail_pass) || !nzchar(mail_user)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Not configured",
        text = "Email sending is not configured on this server. Please report the issue at https://github.com/WangJin93/TFTF/issues instead.",
        type = "warning"
      )
      return(NULL)
    }

    if (!requireNamespace("mailR", quietly = TRUE)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Not available",
        text = "The 'mailR' package is required to send emails but is not installed.",
        type = "error"
      )
      return(NULL)
    }

    mail_to <- Sys.getenv("TFTF_MAIL_TO", unset = getOption("TFTF.mail_to", mail_user))
    ok <- tryCatch({
      mailR::send.mail(from = mail_user,
                       to = mail_to,
                       subject = paste0(input$issue_title, " & ", input$issue_labels),
                       body = input$issue_description,
                       smtp = list(host.name = Sys.getenv("TFTF_MAIL_SMTP", "smtp.qq.com"),
                                   port = as.integer(Sys.getenv("TFTF_MAIL_PORT", "465")),
                                   user.name = mail_user,
                                   passwd = mail_pass,
                                   ssl = TRUE),
                       authenticate = TRUE,
                       send = TRUE)
      TRUE
    }, error = function(e) {
      message("Failed to send feedback email: ", conditionMessage(e))
      FALSE
    })

    # show confirmation
    if (ok) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Issue Submitted!",
        text = "Thank you for your feedback! Your issue has been sent. One of the maintainers will review your issue and contact you for additional details.",
        type = "success"
      )
    } else {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Sending failed",
        text = "Sorry, the email could not be sent. Please report the issue at https://github.com/WangJin93/TFTF/issues instead.",
        type = "error"
      )
    }
    
    # reset key inputs
    updateTextInput(
      session,
      "issue_title",
      value = ""
    )
    
    updateTextAreaInput(
      session,
      "issue_description",
      value = ""
    )
    
    shinyWidgets::updatePickerInput(
      session,
      "issue_labels",
      selected = character(0)
    )
  })

 
}
    
## To be copied in the UI
# mod_feedback_ui("feedback_ui_1")
    
## To be copied in the server
# callModule(mod_feedback_server, "feedback_ui_1")
 
