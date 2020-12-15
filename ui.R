library(shiny)
library(shinyWidgets)
library(shinythemes)

shinyUI(fluidPage(
    theme = shinytheme(theme = "slate"),
    #shinythemes::themeSelector(),
    # setBackgroundColor(
    #     color = "ghostwhite",
    #     gradient = c("linear"),
    #     direction = c("top"),
    #     shinydashboard = TRUE),
    titlePanel("endoGenes - Reference Genes Analysis"),
    sidebarLayout(
        sidebarPanel(
            h5(helpText("The ct file must have a column containing the sample names (as first column), one column per gene, then a groups column (containing the treatment factors) (optional)")),
            uiOutput("example_ct"),
            fileInput("file1", "Upload the Endogenous ct file"),
            radioButtons("group",
                         "Is there a group on the dataset?",

                         choiceNames = list("Yes", "No"),
                         choiceValues = list(TRUE, FALSE),
                         selected = TRUE, inline = TRUE),
            h5(helpText("Is there a header on the endogenous dataset?")),
            checkboxInput(inputId = 'header', label = 'Yes', value = TRUE),
            radioButtons(inputId = 'sep', label = 'Which is the separator?', choices = c(Comma = ',', Semicolon = ';', Tab = '\t', Space = ''), selected = '\t', inline = TRUE),
            hr(),
            hr(),
            h5(helpText("The efficiency file must have a column containing the gene names and one column containing the efficiency values")),
            br(),
            uiOutput("example_eff"),
            br(),
            fileInput("file2", "Upload the efficiency file"),
            h5(helpText("Is there a header on the efficiency dataset?")),
            checkboxInput(inputId = 'header2', label = 'Yes', value = TRUE),
            radioButtons(inputId = 'sep2', label = 'Which is the separator?', choices = c(Comma = ',', Semicolon = ';', Tab = '\t', Space = ''), selected = '\t', inline = TRUE),
            tags$hr(),
            #br(),
            width = 3, 
            hr(),
            hr(),
            hr(),
            uiOutput("links"),
            #h5(helpText("This shinyApp uses the following packages: RankAggreg, ctrlGene, SLqPCR and, a function provided by https://moma.dk/normfinder-software"))
        ),
        mainPanel(
            uiOutput("tb")
            #plotOutput("plt1"),
            
            #plotOutput("plt2"),
            #downloadButton('down', 'Download the Plot')
            
        )
        
    )
))




