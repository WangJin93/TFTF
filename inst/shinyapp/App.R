library(TFTF)
if(!require("shiny")) install.packages("shiny",update = F,ask = F)
if(!require("waiter")) install.packages("waiter",update = F,ask = F)
if(!require("ggplot2")) install.packages("ggplot2",update = F,ask = F)
if(!require("ggpubr")) install.packages("ggpubr",update = F,ask = F)
if(!require("dplyr")) install.packages("dplyr",update = F,ask = F)
if(!require("DT")) install.packages("DT",update = F,ask = F)
if(!require("httr")) install.packages("httr",update = F,ask = F)
if(!require("stringr")) install.packages("stringr",update = F,ask = F)
if(!require("jsonlite")) install.packages("jsonlite",update = F,ask = F)
if(!require("VennDiagram")) install.packages("VennDiagram",update = F,ask = F)
if(!require("shinyWidgets")) install.packages("shinyWidgets",update = F,ask = F)
if(!require("UCSCXenaShiny")) BiocManager::install("UCSCXenaShiny",update = F,ask = F)
if(!require("psych")) install.packages("psych",update = F,ask = F)
if(!require("igraph")) install.packages("igraph",update = F,ask = F)
if(!require("ggraph")) install.packages("ggraph",update = F,ask = F)
if(!require("tidygraph")) install.packages("tidygraph",update = F,ask = F)
if(!require("bs4Dash")) install.packages("bs4Dash",update = F,ask = F)
if(!require("plotrix")) install.packages("plotrix",update = F,ask = F)

source("apps/Dashboard.R")
source("apps/tf_target.R")
source("apps/tf_target_net.R")
source("apps/target_tf.R")
source("apps/tcga_genecor.R")
source("apps/tf_list.R")
source("apps/feedback.R")
# Define UI for application that draws a histogram
ui <- dashboardPage(
  dashboardHeader(
    title = h3("TF-Target Finder", style = 'font-size:28px;color:#ffffff;font-weight: bold;font-family: "Georgia"'),
    skin = "info"),

  dashboardSidebar(status = "info",
    sidebarMenu(
      menuItem("Welcome",icon = icon("info-circle") , tabName = "Welcome"),
      menuItem("TF-->Targets",icon = icon('palette') , tabName ="TF"),
      menuItem("Target-->TFs",icon = icon("palette") , tabName ="Target"),
      menuItem("Pan-tissue cor",icon = icon('deezer') , tabName ="tcga"),
      menuItem("TF-target net",icon = icon('glyphicon glyphicon-fullscreen',lib="glyphicon") , tabName ="Net"),
      menuItem("TF list",icon = icon('th') , tabName ="tf_list"),
      menuItem("Feedback",icon = icon('envelope') , tabName ="feedback")
    )
  ),

  dashboardBody(
    tabItems(
      tabItem("Welcome",ui.modules_dash("Welcome")),
      tabItem("TF",ui.modules_tf("TF")),
      tabItem("Net",ui.modules_net("Net")),
      tabItem("Target",ui.modules_target("Target")),
      tabItem("tcga",ui.modules_tcga("tcga")),
      tabItem("tf_list",ui.modules_tf_list("tf_list")),
      tabItem("feedback",mod_feedback_ui("feedback"))
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  callModule(server.modules_dash, "Welcome")
  callModule(server.modules_tf, "TF")
  callModule(server.modules_net, "Net")
  callModule(server.modules_target, "Target")
  callModule(mod_feedback_server, "feedback")
  callModule(server.modules_tcga, "tcga")
  callModule(server.modules_tf_list, "tf_list")

}

# Run the application
shinyApp(ui = ui, server = server)
