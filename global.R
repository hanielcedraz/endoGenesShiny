library(pacman)

p_load(
    shiny,
    ggplot2,
    RankAggreg,
    ctrlGene,
    reshape2,
    RColorBrewer,
    data.table,
    DT,
    dplyr,
    shinyBS,
    SLqPCR,
    bs4Dash,
    plotly,
    shinycustomloader
    )


source("RnormfinderFunction.txt")



# 
# getCtvalues <- function(inFile1) {
#     #inFile1 <- input$file1
#     if (is.null(inFile1)) {return(NULL)}
#     dataset <- fread(inFile1$datapath, header = input$header, sep = input$sep, data.table = FALSE)
#     colnames(dataset) <- c("SampleName", colnames(dataset[,2:ncol(dataset)]))
#     
#     return(as_tibble(dataset))
# }
# 
# efficiency <- reactive({
#     inFile2 <- input$file2
#     if (is.null(inFile2)) {return(NULL)}
#     fread(inFile2$datapath, header = input$header2, sep = input$sep2, col.names = c('GENE', 'EFFICIENCY'))
# })