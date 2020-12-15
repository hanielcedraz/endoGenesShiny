library(shiny)
library(shinydashboard)
library(ggplot2)
library(RankAggreg)
library(ctrlGene)
library(reshape2)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(shinyBS)
library(SLqPCR)

#getwd()
#do.call(source("./RnormfinderFunction.txt"), )
source("RnormfinderFunction.txt")

#setwd("~/")
# Define server logic required to do the back-end ######################
shinyServer(
  function(input, output) {

    
## defining the reactive functions to sign the basic variables such as dataset and basic plots        
    ctvalues <- reactive({
        inFile1 <- input$file1
        if (is.null(inFile1)) {return(NULL)}
        dataset <- fread(inFile1$datapath, header = input$header, sep = input$sep, data.table = FALSE)
        colnames(dataset) <- c("SampleName", colnames(dataset[,2:ncol(dataset)]))
        
        return(as_tibble(dataset))
    })
    
    efficiency <- reactive({
        inFile2 <- input$file2
        if (is.null(inFile2)) {return(NULL)}
        fread(inFile2$datapath, header = input$header2, sep = input$sep2, col.names = c('GENE', 'EFFICIENCY'))
    })
    
    # output$filedf <- renderTable({
    #     if (is.null(ctvalues())) {return()}
    #     input$file
    # })
    
    # output$sum <- renderTable({
    #     if (is.null(ctvalues())) {return()}
    #     summary(ctvalues())
    # })
    
    output$table <- renderTable({
      if (!is.null(input$file1$datapath)) {
        if (input$head == TRUE) {
          return(ctvalues())
        } 
        else {
          return(head(ctvalues()))
        }
      } else {
        return(paste("Please, upload the endogenous Ct file"))
      } 
    }) 
    

    output$effi <- renderTable({
        if (!is.null(input$file2$datapath)) {
          return(efficiency())
        } else {
          return(paste("Please, upload the Efficiency file"))
        }
    })
    
    reshapeCt <- reactive({
        ctvaluenew <- ctvalues()[,2:ncol(ctvalues())]
        #dat.m <- reshape2::melt(ctvaluenew, measure.vars = colnames(ctvaluenew[-ncol(ctvaluenew)]))
        reshape <- reshape2::melt(ctvaluenew, measure.vars = colnames(ctvaluenew[-ncol(ctvaluenew)]))
        colnames(reshape) <- (c('groups', 'Gene', 'ctValue'))
        return(reshape)
    })
    
    output$rephase <- renderTable({
        if (is.null(reshapeCt())) {return()}
        reshapeCt()
    })
    
    
    boxplot1 <- reactive({
        ggplot(reshapeCt()) +
          geom_boxplot(aes(x = groups, y = ctValue, color = Gene), 
                       outlier.shape = NA, fill = "transparent") +
          labs(x = "Genes", y = "Ct number") +
          scale_fill_brewer(palette = "RdBu") + 
          theme_minimal()

        
    })
    
    boxplot2 <- reactive({
        ggplot(reshapeCt()) +
          geom_boxplot(aes(x = Gene, y = ctValue, color = groups), outlier.shape = NA) +
          labs(x = "Genes/Groups", y = "Ct number") +
          scale_fill_brewer(palette = "RdBu") + 
          theme_minimal()
    })
    
    output$boxplots <- renderPlot({
      if (!is.null(input$file1$datapath)) {
        if (input$group == TRUE) {
          if (input$plot == "gene") {
            boxplot1()
          }
          else if (input$plot == "group") {
            boxplot2()
        }
        } else if (input$group == FALSE) {
          if (input$plot == "gene") {
            boxplot1()
          }
          else if (input$plot == "group") {
            return(paste("Theres no group on the dataset"))
          }
          }
      } else {
        return(paste("Please, upload the endogenous Ct file and the efficieny file"))
        }
    })
    
    output$downloadBoxPlots <- downloadHandler(
      filename = function() {
        if (input$plot == "gene") {
          paste("boxplotByGene", input$type, sep = '.')
        }
        else if (input$plot == "group") {
          paste("boxplotByGroup", input$type, sep = '.')
        }
      },
      content = function(file) {
        if (input$plot == "gene") {
          ggsave(file, plot = boxplot1(), device = input$type)
        }
        else if (input$plot == "group") {
          ggsave(file, plot = boxplot2(), device = input$type)
        }
      }
    )
    
    output$downloadbp2 <- downloadHandler(
        filename = function() {
            paste("boxplot2", input$type2, sep = '.')
        },
        
        content = function(file) {
            ggsave(file, plot = boxplot2(), device = input$type2)
        }
    )
    
############# End of the first block
#########################################################################        
            
############ Normalizing the dataset according to vandesopelen
    
    
    new <- reactive({
      if (input$group == TRUE) {
        news <- ctvalues()[,-length(ctvalues())]
        return(news)
      }
      else if (input$group == FALSE) {
        news <- ctvalues()
        return(news)
      }
      })
    
    qvalue <- reactive({
      inFile1 <- input$file1
      colValues <- read.table(inFile1$datapath, header = input$header, sep = input$sep, row.names = 1)
      if (input$group == TRUE) {
        colvalue_new <- colValues[,-length(colValues)]
        #return(new)
      } else if (input$group == FALSE) {
        colvalue_new <- colValues
        #return(new)
      } #
      
      
      
      
      efficiency()  # Efficiency values
      rel_values <- function(dados, efi_list){
        res <- list()
        for (i in 1:dim(dados)[2]) {
          res[[i]] <- relQuantPCR(dados[,i], E = efi_list[i], na.rm = FALSE)
          
        }
        res_matrix <- do.call('cbind',res);
      }
      
      qvalue <- rel_values(colvalue_new, efficiency()$EFFICIENCY) 
      ##nomes_col <- colnames(c("SampleName", colnames(new)))
      qvalue2 <- data.table::data.table(SampleName = ctvalues()[1], qvalue)
      colnames(qvalue2) <- c("SampleName", colnames(colvalue_new[1:ncol(colvalue_new)]))
      return(qvalue2)
    })
    
   
    output$tables <- renderTable({
      if (input$dataset == "1") {
        if (!is.null(input$file1$datapath) & !is.null(input$file2$datapath)) {
          if (input$head == TRUE) {
            return(ctvalues())
          } 
          else {
            return(head(ctvalues()))
          }
        } 
        else {
          return(paste("Please, upload the endogenous Ct file and the efficieny file"))
        }} 
        else if (input$dataset == "2") {
          if (!is.null(input$file1$datapath) & !is.null(input$file2$datapath)) {
          if (input$head == TRUE) {
            return(qvalue())
          } 
          else {
            return(head(qvalue()))
          }
        }
          else {
            return(paste("Please, upload the endogenous Ct file and the efficieny file"))
          }
        }
      else if (input$dataset == "3") {
        if (!is.null(input$file1$datapath) & !is.null(input$file2$datapath)) {
        if (input$head == TRUE) {
          return(efficiency())
        } 
        else {
          return(head(efficiency()))
        }
      }
        else {
          return(paste("Please, upload the endogenous Ct file and the efficieny file"))
        }
      }
    })
    
    output$downloadNormData <- downloadHandler(
        filename = function() {
            paste("NormData", "csv", sep = ".")
        },
        
        content = function(file) {
            fwrite(qvalue(), file, quote = FALSE, sep = "\t")
        }
    )
    
####################### REFERENCE GENES ANALISYS ######################
    # From here the code is for analisys
    # 1. Genorm
    # 2. NormFinder
    # 3. Bestkeeper
    # 4. Final Ranking
    
#################### STARTING GeNorm analysis ####### 
    
#Using SLqPCR package
## Selection of reference/housekeeping genes

    
    rankingGenorm <- reactive({
      if (input$group == TRUE) {
        ctvalue_new <- ctvalues()[,-length(ctvalues())]
        #return(ctvalue_new)
    } else if (input$group == FALSE) {
        ctvalue_new <- ctvalues()
        #return(ctvalue_new)
    }
      stabGenes <- geNorm(ctvalue_new[,2:ncol(ctvalue_new)], genes = data.frame(Genes = character(0), Avg.M = numeric(0)), ctVal = TRUE)
      stabTabGenes <- data.table(Rank = nrow(stabGenes):1, stabGenes)
      stabTabGenesOrdered <- stabTabGenes[order(stabTabGenes$Avg.M),]
      return(stabTabGenesOrdered)
    })
    
    pairwiseGenorm <- reactive({
      #ctvalue_new <- ctvalues()[,-length(ctvalues())]
      pairwise <- pairwiseV(qvalue()[,2:ncol(qvalue())], ctVal = FALSE)
      colnames(pairwise) <- c("PairWise", "Values")
      return(pairwise)
    })
    
    
##### rendering tab
    
   output$geNormtable <- renderTable({
       if (is.null(rankingGenorm())) {return()}
     if (input$genormtable == "1") {
       rankingGenorm()
     } else if (input$genormtable == "2") {
       pairwiseGenorm()
     }
     
    })

   output$downloadGeNormranking <- downloadHandler(
       filename = function() {
         if (input$genormtable == "1") {
           paste("geNormRanking", "csv", sep = ".")
         } else if (input$genormtable == "2") {
           paste("geNormPairWise", "csv", sep = ".")
         }
       },
       content = function(file) {
         if (input$genormtable == "1") {
           fwrite(rankingGenorm(), file, quote = FALSE, sep = "\t")
         } else if (input$genormtable == "2") {
           fwrite(pairwiseGenorm(), file, quote = FALSE, sep = "\t")
         }
       }
   )
   
   
    ## Ranking Stability Gene
    
   gnormPlotRanking <- reactive({
     ggplot(rankingGenorm(), 
            aes(x = Rank, y = Avg.M, color = Genes)) +
       geom_line(color = input$color) +
       geom_point() +
       lims(x = rankingGenorm()$Gene, y = c(0, max(rankingGenorm()$Avg.M))) +
       ylab("Average expression stability M") +
       xlab("<==== Most Stable Gene   Least Stable Gene ====>") +
       ggtitle("Gene stability measure by ctrlGene Package (geNorm)") +
       theme(
         plot.title = element_text(size = 14, hjust = 0.5)) +
       scale_fill_brewer(palette = "RdBu") + 
       theme_minimal()
   })
   

   pairwisePlotGenorm <- function(){
     #ctvalue_new <- ctvalues()[,-length(ctvalues())]
     pairwise <- pairwiseV(qvalue()[,2:ncol(qvalue())], ctVal = FALSE)
     colnames(pairwise) <- c("PairWise", "Values")
     mypalette <- brewer.pal(8, "Spectral")
     barplot(cbind(pairwise$Values), beside = TRUE, col = mypalette, space = c(0, 0.5), ylim = c(0, 0.25), xlab = c("PairWise Variations"))
     axis(1, at = 1:length(pairwise$Values), labels = as.character(pairwise$PairWise))
     abline(h = seq(0, max(pairwise$Values), by = 0.2), lty = 2, col = "grey")
     abline(h = 0.15, lty = 1, col = input$color)
     box()
   }
   
   
   output$plotGenormShow <- renderPlot({
     if (input$genormplot == "1") {
       gnormPlotRanking() 
     } else if (input$genormplot == "2") {
       pairwisePlotGenorm()
     }
      
     #gnormPlotRanking()
   })
   
   
   output$downloadgenormplots <- downloadHandler(
       filename = function() {
         if (input$genormplot == "1") {
           paste("genormStabilityPlot", input$type3, sep = '.')
         } else if (input$genormplot == "2") {
           paste("genormPairWisePlot", input$type3, sep = '.')
         }
       },
       content = function(file) {
         if (input$genormplot == "1") {
         ggsave(file, plot = gnormPlotRanking(), device = input$type3)
         } else if (input$genormplot == "2") {
           if (input$type3 == "png") 
             png(file)
            else 
             pdf(file)
             pairwisePlotGenorm()
             dev.off()
           
       }
         })
    

   output$downloadpairwise <- downloadHandler(
     filename = function() {
       paste("pairWisePlot", input$type4, sep = ".")
     },
     content = function(file) {
       if (input$type4 == 'png') 
         png(file)
       else
         pdf(file)
       pairwisePlotGenorm()
       dev.off()
     }
   )
   
################# NORMFINDER #################
   
###### normfind script block
   
   
   
   normTrans <- reactive({
     inFile1 <- input$file1
     if (is.null(inFile1)) {return(NULL)}
     ctvalues2 <- read.table(inFile1$datapath, header = input$header, sep = input$sep, row.names = 1)
     ctvalueTrans <- t(ctvalues2)

       return(ctvalueTrans)
     
   })
   
   output$transtable <- renderTable({
       normTrans()
   },  rownames = TRUE)
   
   output$downloadTransNorm <- downloadHandler(
       filename = function() {
           paste("dataTransNormFinder", "txt", sep = '.')
       },
       content = function(file) {
           write.table(normTrans(), file)
       }
   )
   
   normFinderFunction <- reactive({
        inFile3 <- input$file3
       if (is.null(inFile3)) {return(NULL)}
       # ctvalues2 <- read.table(inFile3$datapath, header = input$header, sep = input$sep)
       
       
       resulttotal <- Normfinder(inFile3$datapath, Groups = input$group, ctVal = TRUE)
       return(resulttotal)
       
       
   }) # script from normfinder algorithm
   
  
   stabilityNofmFinder <- reactive({
     if (input$group == TRUE) {
       stabilityG <- data.table(rank = 1:nrow(normFinderFunction()$Ordered), normFinderFunction()$Ordered, keep.rownames = TRUE)
       colnames(stabilityG) <- c("Rank", "Gene", "GroupDif", "GroupSD", "Stability")
       return(stabilityG)
     } else if (input$group == FALSE) {
       stabilityNG <- data.table(rank = 1:nrow(normFinderFunction()$Ordered), Gene = rownames(normFinderFunction()$Ordered), normFinderFunction()$Ordered$GroupSD)
       colnames(stabilityNG) <- c("Rank", "Gene", "Stability")
       return(stabilityNG)
       }
   })
   
   # 
   # stabilityNofmFinderNG <- reactive({
   #   stabilityNG <- data.table(rank = 1:nrow(normFinderFunction()$Ordered), Gene = rownames(normFinderFunction()$Ordered), normFinderFunction()$Ordered$GroupSD)
   #   colnames(stabilityNG) <- c("Rank", "Gene", "Stability")
   #   return(stabilityNG)
   # })
   
   
   output$normFinderTable <- renderTable({
     if (!is.null(input$file3$datapath)) {
       if (input$group == TRUE) {
         if (input$normfindertable == "1") {
            stabilityNofmFinder()
          } else if (input$normfindertable == "2") {
             normFinderFunction()$PairOfGenes
        }
          } else if (input$group == FALSE) {
            if (input$normfindertable == "2") {
              return(paste("This table will be generated only if there is a groups column in dataset"))
            } else {
              stabilityNofmFinder()
            }
            
        #return(paste("Groups igual FALSE"))
          } 
       } else {
       return(paste("Please, download the file 'dataTransNormFinder.txt' on the previous tab and upload it again"))
     }
   })


   output$dowRankNormFinder <- downloadHandler(
     filename = function() {
       paste("rankNormfinder", "csv", sep = ".")
     },
     content = function(file) {
       fwrite(stabilityNofmFinder(), file, quote = FALSE, sep = "\t")
     }
   )
   
    normfinderPlot <- reactive({
      ggplot(stabilityNofmFinder(), aes(x = Rank, y = Stability, color = Gene)) +
        geom_line(color = input$color2) +
        geom_point() +
        lims(x = stabilityNofmFinder()$Gene, y = c(0, max(stabilityNofmFinder()$Stability))) +
        ylab("Average expression stability (S)") +
        xlab("<==== Most Stable Gene   Least Stable Gene ====>") +
        ggtitle("Gene Stability Measure by NormFinder") +
        theme(
          plot.title = element_text(size = 14, hjust = 0.5)) +
        scale_fill_brewer(palette = "RdBu") + 
        theme_minimal()
        
    })
   
   output$normFinderPlot <- renderPlot({
     if (!is.null(input$file3$datapath)) {
       normfinderPlot()
     } else {
       return(paste("Please, download the file 'dataTransNormFinder.txt' on the previous tab and upload it again"))
     }
   })
   
   output$downnormfinderplot <- downloadHandler(
     filename = function() {
       paste("normFinderStabilityPlot", input$type5, sep = ".")
     },
     content = function(file) {
       ggsave(file, normfinderPlot(), device = input$type5)
     }
   )
 
   
   
########################### BESTKEEPER ############################
   
   bestkeepanaly <- reactive({
     inFile1 <- input$file1
     if (is.null(inFile1)) {return(NULL)}
     ctvalue <- read.table(inFile1$datapath, header = input$header, sep = input$sep, row.names = 1)
     if (input$group == TRUE) {
     bestkeeper_results = bestKeeper(ctvalue[,1:ncol(ctvalue) - 1], ctVal = TRUE)
     return(bestkeeper_results)
     } else if (input$group == FALSE) {
       bestkeeper_results = bestKeeper(ctvalue[,1:ncol(ctvalue)], ctVal = TRUE)
       return(bestkeeper_results)
     }
   })
   
   cpstatbestkeeper <- reactive({
     cpStatBestkeeper <- data.table(bestkeepanaly()$CP.statistics, keep.rownames = TRUE)
     return(cpStatBestkeeper)
   })
   
   pair.Wise.corbestkeeper <- reactive({
     pairwisebestkeeper <- data.table(bestkeepanaly()$pair.Wise.cor, keep.rownames = TRUE)
     return(pairwisebestkeeper)
   })
   
   HKG.vs.BestKeeper <- reactive({
     hkgvsbestkeeper <- data.table(bestkeepanaly()$HKG.vs.BestKeeper, keep.rownames = TRUE)
     return(hkgvsbestkeeper)
   })
   
   finalbestrank <- reactive({
     bestkeeper_genes <- t(sort(bestkeepanaly()$CP.statistics[6,]));
     best_genes <- t(bestkeeper_genes);
     best_genes_values <- data.table(Rank = 1:nrow(best_genes), rownames(best_genes), best_genes)
     colnames(best_genes_values) <- c('Rank', 'Genes', 'SD_Value')
     return(best_genes_values)
   })
   
    output$bestkeepertables <- renderTable({
      if (input$table == "1") {
        cpstatbestkeeper()
        #as.data.frame(bestkeepanaly()$CP.statistics)
        }
      else if (input$table == "2") {
        pair.Wise.corbestkeeper()
        #bestkeepanaly()$pair.Wise.cor
        }
      else if (input$table == "3") {
        HKG.vs.BestKeeper()
        #bestkeepanaly()$HKG.vs.BestKeeper
      }
      else if (input$table == "4") {
        finalbestrank()
}
})
          
   
    output$dowbestkeeper <- downloadHandler(
      filename = function() {
        if (input$table == "1") {
          paste("CPstatisticsBestkeeper", "csv", sep = ".")
        }
        else if (input$table == "2") {
          paste("pairWisecorBestkeeper", "csv", sep = ".")
        }
        else if (input$table == "3") {
          paste("HKG.vs.BestKeeperNestkeeper", "csv", sep = ".")
        }
        else if (input$table == "4") {
        paste("finalRankingBestkeeper", "csv", sep = ".")
        }
      }
        ,
      content = function(file) {
        if (input$table == "1") {
          fwrite(cpstatbestkeeper(), file, quote = FALSE)
        }
        else if (input$table == "2") {
          fwrite(pair.Wise.corbestkeeper(), file, quote = FALSE, sep = "\t")
        }
        else if (input$table == "3") {
          fwrite(HKG.vs.BestKeeper(), file, quote = FALSE, sep = "\t")
        }
        else if (input$table == "4") {
          fwrite(finalbestrank(), file, quote = FALSE, sep = "\t")
        }
      }
    )
   
    
    plotreacbestkeepr <- reactive({
      ggplot(finalbestrank(), 
             aes(x = Rank, y = SD_Value, color = Genes)) +
        geom_line(color = input$color3) +
        geom_point() +
        lims(x = finalbestrank()$Genes, y = c(0, max(finalbestrank()$SD_Value))) +
        ylab("Average expression stability (Standard Deviation - SD)") +
        xlab("<==== Most Stable Gene   Least Stable Gene ====>") +
        ggtitle("Gene Stability Measure by Bestkeeper") +
        theme(
          plot.title = element_text(size = 14, hjust = 0.5)) +
        scale_fill_brewer(palette = "RdBu") + 
        theme_minimal()
    })
    

    output$plotbestkeeper <- renderPlot({
      plotreacbestkeepr()
    })
    
    output$dowFinalRankBestkeeper <- downloadHandler(
      filename = function() {
        paste("bestkeeperPlot", input$type6, sep = ".")
      },
      content = function(file) {
        ggsave(file, plotreacbestkeepr(), device = input$type6)
      }
    )
    
    
    
    
    CESP <- reactive({
      genes <- strsplit(rankingGenorm()$Genes[1], split = "-")
      # gene1 <- substr(basename(paste0(rankingGenorm()$Genes[1])), 1, nchar(basename(paste0(rankingGenorm()$Genes[1]))) - 7)
      # gene2 <- substr(basename(paste0(rankingGenorm()$Genes[1])), 8, nchar(basename(paste0(rankingGenorm()$Genes[1]))) - 0)
      
      ratools <- t(data.frame(NormFinder = stabilityNofmFinder()$Gene, GeNorm = c(genes[[1]][1], genes[[1]][2], rankingGenorm()$Genes[2:length(rankingGenorm()$Genes)]), Bestkeeper = finalbestrank()$Genes))
      colnames(ratools) <- paste(1:ncol(ratools))
      
      
      ra_weight_ratools <- t(data.frame(NormFinder = stabilityNofmFinder()$Stability, GeNorm = c(rankingGenorm()$Avg.M[1], rankingGenorm()$Avg.M), Bestkeeper = finalbestrank()$SD_Value))
      colnames(ra_weight_ratools) <- paste(1:ncol(ratools))
      
      
      k <- ncol(ratools)
      
      #method = opt$method
      #distance = opt$distance
      #maxIter = opt$iteraction
      
      
      if (k > 10) {
        CESP <- RankAggreg(ratools, k = ncol(ratools), ra_weight_ratools, method = input$method, distance = input$distance, weight = .25, rho = .1, maxIter = 1000, verbose = TRUE)
      }else {
        CESP <- BruteAggreg(ratools, k = ncol(ratools), ra_weight_ratools, distance = input$distance)
      }
      
      return(CESP)
    })
    
    rankfinal <- reactive({
      final_res <- data.frame(Rank = paste(1:length(CESP()$top.list)), Gene = CESP()$top.list)
      return(final_res)
    })
    
    output$finalranking <- renderTable({
      rankfinal()
    })
    
    
    output$dowfinalrank <- downloadHandler(
      filename = function() {
        paste("finalRanking", "csv", sep = ".")
      },
      content = function(file) {
        fwrite(rankfinal(), file, quote = FALSE, sep = "\t")
      }
    )
    
    output$plotfinal <- renderPlot({
      plot(CESP())
    })
    
    output$dowfinalPlot <- downloadHandler(
      filename = function() {
        paste("finalRankingPlot", input$type7, sep = ".")
      },
      content = function(file) {
        if (input$type7 == 'png') 
          png(file)
        else 
          pdf(file)
          plot(CESP())
          dev.off()
        
      }
    )

 
    url_ct <- a("example txt file with ct data.", href = "https://github.com/hanielcedraz/endoGenes/blob/master/endogenous_ct.txt")
    
    url_eff <- a("example txt file with efficiency data.", href = "https://github.com/hanielcedraz/endoGenes/blob/master/efficiencies_list.txt")
    
    output$example_ct <- renderUI({
      tagList("Follow this link to download an", url_ct)
    })
    
    output$example_eff <- renderUI({
      tagList("Follow this link to download an", url_eff)
    })

    RankAggreg <- a("RankAggreg, ", href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-62")
    ctrlGene <- a("ctrlGene, ", href = "https://cran.r-project.org/web/packages/ctrlGene/index.html")
    SLqPCR <- a("SLqPCR, ", href = "https://www.bioconductor.org/packages/release/bioc/html/SLqPCR.html")
    normfinder <- a("normfinder", href = "https://moma.dk/normfinder-software")
    
    output$links <- renderUI({
      tagList("This shinyApp uses the following packages: ", RankAggreg, ctrlGene, SLqPCR, "and, a function for", normfinder, "provided by MOMA -Department of Molecular Medicine")
    })
    
    
    output$tb <- renderUI({
        dashboardSidebar()
        dashboardBody(
          # Boxes need to be put in a row (or column)
        tabsetPanel(
            tabPanel("DataSets",
                     #h5(helpText("Is there a group on the dataset?")),
                     # radioButtons("group", 
                     #              "Is there a group on the dataset?",
                     #              
                     #              choiceNames = list("Yes", "No"),
                     #              choiceValues = list(TRUE, FALSE),
                     #              selected = TRUE, inline = TRUE),
                     # checkboxInput(inputId = 'group', label = 'Yes',
                     #               value = TRUE),
               selectInput("dataset", "Select the Table",
                           choices = c("Endogenous Dataset" = "1",
                                       "Normalized Data" = "2",
                                       "Efficiency File" = "3"
                           )),
               tableOutput("tables"),
               tipify(checkboxInput(inputId = 'head', 
                                    label = "Expand", 
                                    value = FALSE),
                      "Click here to show the entire dataset",
                      placement = "left", trigger = "hover"),
               downloadButton("downloadNormData", 
                              "Download the Norm data")),
              # tableOutput("effi")),
      #tabPanel("Summary", tableOutput("sum")),
      #tabPanel("Efficiency", tableOutput("effi")),
      #tabPanel("reshape", tableOutput("rephase")),
      tabPanel("Boxplot ct Values",
               selectInput("plot", "Select the Plot",
                           choices = c("Boxplot by Gene" =
                                         "gene",
                                       "Boxplot by Group" =
                                         "group")),
               plotOutput("boxplots"),
               downloadButton("downloadBoxPlots",
                              "Download the Plot"),
               radioButtons("type", "Choose the file type",
                            choices = list('png', 'pdf', "tiff"),
                            selected = 'png')),
      tabPanel("Ct values transposed", 
               h5(helpText("This dataset will be used
                         in the next analysis. Please, 
                           download it")),
               tableOutput("transtable"),
               downloadButton("downloadTransNorm", "Download
                              the transposed dataset"),
               fileInput("file3","Please, download the above
                           dataset and upload it again
                         in this box. This dataset will be used
                         in the next analysis")),
      tabPanel("NormFinder Results",
               box(selectInput("normfindertable", "Select the Table",
                               choices = c("NormFinder Ranking" = "1",
                                           "Pair of Genes" = "2"
                               )),
                 tableOutput("normFinderTable"),
                   # h5(helpText("Is there a group on the dataset?")),
                   # checkboxInput(inputId = 'group', label = 'Yes',
                   #               value = TRUE),
                   downloadButton("dowRankNormFinder", "Download
                              the file")),
               box(plotOutput("normFinderPlot"),
                   downloadButton("downnormfinderplot", "Download
                              the file"),
                   textInput("color2", "Choose the line color", 
                             value = "black", placeholder = "Type a valid color"),
                   radioButtons("type5", "Choose the file type",
                                choices = list('png', 'pdf', "tiff"),
                                selected = 'png'), width = 12)),
      tabPanel("geNorm Results",
         box(selectInput("genormtable", "Select the Table",
                     choices = c("Genorm Ranking" = "1",
                                 "PairWise Correlation" = "2"
                                 )),
         tableOutput("geNormtable"),
         downloadButton("downloadGeNormranking",
                        "Download the geNorm Ranking")),
         box(selectInput("genormplot", "Select the Plot",
                         choices = c("Genorm Ranking" = "1",
                                     "PairWise Correlation" = "2"
                         )),
             plotOutput("plotGenormShow"),
             downloadButton("downloadgenormplots", 
                            "Download the Plot"),
             textInput("color", "Choose the line color", 
                       value = "black", placeholder = "Type a valid color"),
             radioButtons("type3", "Choose the file type",
                        choices = list('png', 'pdf', "tiff"),
                        selected = 'png'), width = 12)),
      tabPanel("Bestkeeper Results",
         selectInput("table", "Select the table",
                    choices = c("Statistic Summary" = "1",
                                 "PairWise Correlation" = "2",
                                 "HKG vs. Bestkeeper" = "3",
                                 "Bestkeeper Ranking" = "4")),
         tableOutput("bestkeepertables"),
         downloadButton("dowbestkeeper", "Download the file"),
         plotOutput("plotbestkeeper"),
         textInput("color3", "Choose the line color", 
                   value = "black", placeholder = "Type a valid color"),
         downloadButton("dowFinalRankBestkeeper",
             "Download the file"),
         radioButtons("type6", "Choose the file type",
                    choices = list('png', 'pdf', "tiff"),
                    selected = 'png')),
      tabPanel("Final Ranking", 
         h5(helpText("Please wait until the results to be generate. It could take some time")),
         tableOutput("finalranking"),
         tipify(radioButtons("method", "Choose the method",
                      choices = list("CE", "GA"),
                      selected = "CE", inline = TRUE),
                "Cross Entropy Monte Carlo (CE) or Genetic Algorithm (GA)",
                placement = "left", trigger = "hover"),
         tipify(radioButtons("distance", "Choose the
                             distance",
                      choices = list("Spearman", "Kendall"),
                      selected = "Spearman", inline = TRUE),
                "Distance to be used in the similarity",
                placement = "left", trigger = "hover"),
         downloadButton("dowfinalrank", "Download the file"),
         plotOutput("plotfinal"),
         #h5(helpText("Please wait until the plot to be generate. It could take some time")),
         downloadButton("dowfinalPlot", "Download
                        the file"),
         radioButtons("type7", "Choose the file type",
                      choices = list('png', 'pdf', "tiff"),
                      selected = 'png'))
      # tabPanel("Relative Expression",
      #          #div(style = "display:inline.block",
      #         splitLayout(
      #          fileInput("file4", "Upload the Target Genes ct file"),
      #          fileInput("file5", "Upload the Target Genes efficiency file")),
      #          selectInput("data", "Select the table",
      #                 choices = c("Target Genes Ct" = "1",
      #                             "Target Genes Efficiency" = "2",
      #                             "Target and Reference Genes" = "3")),
      #          tipify(radioButtons("Normethod", "Choose the method",
      #                              choices = list("Vandesompele",
      #                                             "DeltaCt",
      #                                             "DeltaDeltaCt"),
      #                              selected = "DeltaCt", inline = TRUE),
      #                 title = "Choose the method used for normalization",
      #                 placement = "left", trigger = "hover"),
      #          tableOutput("relExp"),
      #          downloadBttn("downRelExp", "Download the table", 
      #                       style = "jelly")
      #          )
        ))
    })
  }
)   
