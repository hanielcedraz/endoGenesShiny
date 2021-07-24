library(shiny)
library(shinyWidgets)
library(shinythemes)
library(fresh)
# library(shinydashboardPlus)

shinyUI(
    ui = dashboardPage(
        title = "endoGenes",
        header = dashboardHeader(
            skin = "light", 
            status = "navy",
            title = dashboardBrand(
                title = "endoGenes - Reference Genes Analysis", 
                color = "navy",
                #image = "stlogo.png",
                opacity = 0.8
                )
            ),
        sidebar = dashboardSidebar(
            expandOnHover = TRUE,
            skin = "light", 
            minified = FALSE, 
            status = "navy",
            sidebarMenu(
                id = "sidebar",
                uiOutput("example_ct"),
                fileInput("file1", "Upload the Endogenous ct file"),
                radioButtons("group",
                             "Is there a group on the dataset?",
                             choiceNames = list("Yes", "No"),
                             choiceValues = list(TRUE, FALSE),
                             selected = TRUE, 
                             inline = TRUE
                             ),
                #h5(helpText("Is there a header on the endogenous dataset?")),
                checkboxInput(
                    inputId = 'header', 
                    label = 'Yes', 
                    value = TRUE
                    ),
                radioButtons(
                    inputId = 'sep', 
                    label = 'What is the separator?', 
                    choices = c("Comma" = ',', 
                                "Semicolon" = ';', 
                                "Tab" = '\t', 
                                "Space" = ''
                                ), 
                    selected = '\t', 
                    inline = TRUE
                    ),
                hr(),
                hr(),
                h5(helpText("The efficiency file must have a column containing the gene names and one column containing the efficiency values")),
                br(),
                uiOutput("example_eff"),
                br(),
                fileInput("file2", "Upload the efficiency file"),
                #h5(helpText("Is there a header on the efficiency dataset?")),
                checkboxInput(
                    inputId = 'header2', 
                    label = 'Yes', 
                    value = TRUE
                    ),
                radioButtons(
                    inputId = 'sep2', 
                    label = 'What is the separator?', 
                    choices = c(
                        "Comma" = ',', 
                        "Semicolon" = ';', 
                        "Tab" = '\t', 
                        "Space" = ''
                        ), 
                    selected = '\t', 
                    inline = TRUE
                    ),
                #tags$hr(),
                #br(),
                width = 3, 
                hr(),
                hr(),
                hr(),
                uiOutput("links")
            )
        ), 
        body = dashboardBody(
            skin = "light",
            status = "danger",
            use_theme(
                create_theme(
                    bs4dash_layout(
                        sidebar_width = "360px"
                    )#,
                    # bs4dash_button(
                    #     default_background_color = "navy",
                    #     default_border_color = "#111"
                    # )
                )
            ),
            useShinyjs(),
            #useShinyalert(),
            fluidPage(
                uiOutput(
                    outputId = "tb"
                    )
                )
            )
        )
    )




