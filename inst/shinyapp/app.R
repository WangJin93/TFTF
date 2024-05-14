#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(waiter)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(DT)
library("httr")
library("stringr")
library("jsonlite")
library(VennDiagram)
library(UCSCXenaShiny)
library(RMySQL)
# library(shinyWidgets)
library(psych)
library(igraph)
library(ggraph)
library(tidygraph)
library(bs4Dash)
source("apps/Dashboard.R")
source("apps/tf_target.R")
source("apps/tf_target_net.R")
source("apps/target_tf.R")
source("apps/tcga_genecor.R")
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
      menuItem("Target-->TFs",icon = icon("th") , tabName ="Target"),
      menuItem("Pan-tissue cor",icon = icon('deezer') , tabName ="tcga"),
      menuItem("TF-target net",icon = icon('glyphicon glyphicon-fullscreen',lib="glyphicon") , tabName ="Net"),
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

}

# Run the application
shinyApp(ui = ui, server = server)
