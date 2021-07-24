
# Define server logic required to do the back-end ######################
shinyServer(function(input, output, session) {


# ## defining the reactive functions to sign the basic variables such as dataset and basic plots
#     ctvalues <- reactive({
#         inFile1 <- input$file1
#         if (is.null(inFile1)) {return(NULL)}
#         dataset <- fread(inFile1$datapath, header = input$header, sep = input$sep, data.table = FALSE)
#         colnames(dataset) <- c("SampleName", colnames(dataset[,2:ncol(dataset)]))
# 
#         return(dataset)
#     })
# 
#     efficiency <- reactive({
#         inFile2 <- input$file2
#         if (is.null(inFile2)) {return(NULL)}
#         eff <- fread(inFile2$datapath, header = input$header2, sep = input$sep2, col.names = c('GENE', 'EFFICIENCY'))
#         return(eff)
#     })
# 
#     # output$filedf <- renderTable({
#     #     if (is.null(ctvalues())) {return()}
#     #     input$file
#     # })
# 
#     # output$sum <- renderTable({
#     #     if (is.null(ctvalues())) {return()}
#     #     summary(ctvalues())
#     # })
# 
#     # output$table <- renderDataTable({
#     #   if (!is.null(input$file1$datapath)) {
#     #     if (input$head == TRUE) {
#     #       return(datatable(ctvalues()))
#     #     }
#     #     else {
#     #       return(head(datatable(ctvalues())))
#     #     }
#     #   } else {
#     #     return(paste("Please, upload the endogenous Ct file"))
#     #   }
#     # })
# 
# 
#     output$effi <- renderTable({
#         if (!is.null(input$file2$datapath)) {
#           return(efficiency())
#         } else {
#           return(paste("Please, upload the Efficiency file"))
#         }
#     })
# 
#     reshapeCt <- reactive({
#         ctvaluenew <- ctvalues()[,2:ncol(ctvalues())]
#         #dat.m <- reshape2::melt(ctvaluenew, measure.vars = colnames(ctvaluenew[-ncol(ctvaluenew)]))
#         if (input$group) {
#           reshape <- reshape2::melt(ctvaluenew, measure.vars = colnames(ctvaluenew[-ncol(ctvaluenew)]))
#           colnames(reshape) <- (c('groups', 'Gene', 'ctValue'))
#           return(reshape)
#         } else {
#           reshape <- reshape2::melt(ctvaluenew, measure.vars = colnames(ctvaluenew))
#           colnames(reshape) <- (c('Gene', 'ctValue'))
#           return(reshape)
#         }
#     })
# 
#     output$rephase <- renderTable({
#         if (is.null(reshapeCt())) {return()}
#         reshapeCt()
#     })
# 
# 
#     boxplot1 <- reactive({
#         ggplot(reshapeCt()) +
#           geom_boxplot(aes(x = groups, y = ctValue, color = Gene),
#                        outlier.shape = NA, fill = "transparent") +
#           labs(x = "Groups", y = "Ct number") +
#           scale_y_continuous(breaks = seq(0, max(reshapeCt()$ctValue), by = 3)) +
#           scale_fill_brewer(palette = "RdBu") +
#           theme_minimal()
# 
# 
#     })
# 
#     boxplotNoGroup <- reactive({
#       ggplot(reshapeCt()) +
#         geom_boxplot(aes(x = Gene, y = ctValue, color = Gene),
#                      outlier.shape = NA, fill = "transparent") +
#         labs(x = "Genes", y = "Ct number") +
#         scale_y_continuous(breaks = seq(0, max(reshapeCt()$ctValue), by = 3)) +
#         scale_fill_brewer(palette = "RdBu") +
#         theme_minimal()
#     })
# 
#     boxplot2 <- reactive({
#         ggplot(reshapeCt()) +
#           geom_boxplot(aes(x = Gene, y = ctValue, color = groups), outlier.shape = NA) +
#           labs(x = "Genes", y = "Ct number") +
#           scale_x_discrete(breaks = seq(0, max(reshapeCt()$ctValue), by = 3)) +
#           scale_fill_brewer(palette = "RdBu") +
#           theme_minimal()
#     })
# 
#     output$boxplots <- renderPlot({
#       if (!is.null(input$file1$datapath)) {
#         if (input$group) {
#             if (input$plot == "gene") {
#               boxplot1()
#             }
#             else if (input$plot == "group") {
#               boxplot2()
#         }
#         } else {
#             if (input$plot == "gene") {
#               boxplotNoGroup()
#             }
#             else if (input$plot == "group") {
#               boxplotNoGroup()
#             }
#           }
#       } else {
#         count.data <- data.frame(
#           class = c("Other reason", "Did not upload the dataset"),
#           n = c(10, 1000)
#         )
# 
#         ggplot(count.data, aes(x = "", y = n, fill = class)) +
#           geom_bar(width = 1, stat = "identity", color = "white") +
#           coord_polar("y", start = 0) +
#           theme_void() +
#           ggtitle(label = "Why is your plot not here?") +
#           theme(plot.title = element_text(color = "black", size = 25, face = "bold.italic", hjust = 0.5, vjust = -5.0), legend.text = element_text(colour = 'black', size = 15), legend.title = element_text(size = 15)) +
#           labs(fill = "Reason")
#         }
#     })
# 
#     output$downloadBoxPlots <- downloadHandler(
#       filename = function() {
#         if (input$plot == "gene") {
#           paste("boxplotByGene", input$type, sep = '.')
#         }
#         else if (input$plot == "group") {
#           paste("boxplotByGroup", input$type, sep = '.')
#         }
#       },
#       content = function(file) {
#         if (input$group) {
#           if (input$plot == "gene") {
#             ggsave(file, plot = boxplot1(), device = input$type)
#           }
#           else if (input$plot == "group") {
#             ggsave(file, plot = boxplot2(), device = input$type)
#           }
#         } else {
#           ggsave(file, plot = boxplotNoGroup(), device = input$type)
#         }
#       }
#     )
# 
#     output$downloadbp2 <- downloadHandler(
#         filename = function() {
#             paste("boxplot2", input$type2, sep = '.')
#         },
# 
#         content = function(file) {
#             ggsave(file, plot = boxplot2(), device = input$type2)
#         }
#     )
# 
# ############# End of the first block
# #########################################################################
# 
# ############ Normalizing the dataset according to vandesopelen
# 
# 
#     new <- reactive({
#       if (input$group == TRUE) {
#         news <- ctvalues()[,-length(ctvalues())]
#         return(news)
#       }
#       else if (input$group == FALSE) {
#         news <- ctvalues()
#         return(news)
#       }
#       })
# 
#     qvalue <- reactive({
#       inFile1 <- input$file1
#       colValues <- read.table(inFile1$datapath, header = input$header, sep = input$sep, row.names = 1)
#       if (input$group == TRUE) {
#         colvalue_new <- colValues[,-length(colValues)]
#         #return(new)
#       } else if (input$group == FALSE) {
#         colvalue_new <- colValues
#         #return(new)
#       } #
# 
# 
# 
# 
#       efficiency()  # Efficiency values
#       rel_values <- function(dados, efi_list){
#         res <- list()
#         for (i in 1:dim(dados)[2]) {
#           res[[i]] <- relQuantPCR(dados[,i], E = efi_list[i], na.rm = FALSE)
# 
#         }
#         res_matrix <- do.call('cbind',res);
#       }
# 
#       qvalue <- rel_values(colvalue_new, efficiency()$EFFICIENCY)
#       ##nomes_col <- colnames(c("SampleName", colnames(new)))
#       qvalue2 <- data.table::data.table(SampleName = ctvalues()[1], qvalue)
#       colnames(qvalue2) <- c("SampleName", colnames(colvalue_new[1:ncol(colvalue_new)]))
#       return(qvalue2)
#     })
# 
# 
#     output$tables <- renderTable({
#       if (input$dataset == "1") {
#         if (!is.null(input$file1$datapath) & !is.null(input$file2$datapath)) {
#           if (input$head == TRUE) {
#             return(ctvalues())
#           }
#           else {
#             return(head(ctvalues()))
#           }
#         }
#         else {
#           return(paste("Please, upload the endogenous Ct file and the efficieny file"))
#         }}
#         else if (input$dataset == "2") {
#           if (!is.null(input$file1$datapath) & !is.null(input$file2$datapath)) {
#           if (input$head == TRUE) {
#             return(qvalue())
#           }
#           else {
#             return(head(qvalue()))
#           }
#         }
#           else {
#             return(paste("Please, upload the endogenous Ct file and the efficieny file"))
#           }
#         }
#       else if (input$dataset == "3") {
#         if (!is.null(input$file1$datapath) & !is.null(input$file2$datapath)) {
#         if (input$head == TRUE) {
#           return(efficiency())
#         }
#         else {
#           return(head(efficiency()))
#         }
#       }
#         else {
#           return(paste("Please, upload the endogenous Ct file and the efficieny file"))
#         }
#       }
#     })
# 
#     output$downloadNormData <- downloadHandler(
#         filename = function() {
#             paste("NormData", "csv", sep = ".")
#         },
# 
#         content = function(file) {
#             fwrite(qvalue(), file, quote = FALSE, sep = "\t")
#         }
#     )
# 
# ####################### REFERENCE GENES ANALISYS ######################
#     # From here the code is for analisys
#     # 1. Genorm
#     # 2. NormFinder
#     # 3. Bestkeeper
#     # 4. Final Ranking
# 
# #################### STARTING GeNorm analysis #######
# 
# #Using SLqPCR package
# ## Selection of reference/housekeeping genes
# 
# 
#     rankingGenorm <- reactive({
#       if (input$group == TRUE) {
#         ctvalue_new <- ctvalues()[,-length(ctvalues())]
#         #return(ctvalue_new)
#     } else if (input$group == FALSE) {
#         ctvalue_new <- ctvalues()
#         #return(ctvalue_new)
#     }
#       stabGenes <- geNorm(ctvalue_new[,2:ncol(ctvalue_new)], genes = data.frame(Genes = character(0), Avg.M = numeric(0)), ctVal = TRUE)
#       stabTabGenes <- data.table(Rank = nrow(stabGenes):1, stabGenes)
#       stabTabGenesOrdered <- stabTabGenes[order(stabTabGenes$Avg.M),]
#       return(stabTabGenesOrdered)
#     })
# 
#     pairwiseGenorm <- reactive({
#       #ctvalue_new <- ctvalues()[,-length(ctvalues())]
#       pairwise <- pairwiseV(qvalue()[,2:ncol(qvalue())], ctVal = FALSE)
#       colnames(pairwise) <- c("PairWise", "Values")
#       return(pairwise)
#     })
# 
# 
# ##### rendering tab
# 
#    output$geNormtable <- renderTable({
#        if (is.null(rankingGenorm())) {return()}
#      if (input$genormtable == "1") {
#        rankingGenorm()
#      } else if (input$genormtable == "2") {
#        pairwiseGenorm()
#      }
# 
#     })
# 
#    output$downloadGeNormranking <- downloadHandler(
#        filename = function() {
#          if (input$genormtable == "1") {
#            paste("geNormRanking", "csv", sep = ".")
#          } else if (input$genormtable == "2") {
#            paste("geNormPairWise", "csv", sep = ".")
#          }
#        },
#        content = function(file) {
#          if (input$genormtable == "1") {
#            fwrite(rankingGenorm(), file, quote = FALSE, sep = "\t")
#          } else if (input$genormtable == "2") {
#            fwrite(pairwiseGenorm(), file, quote = FALSE, sep = "\t")
#          }
#        }
#    )
# 
# 
#     ## Ranking Stability Gene
# 
#    gnormPlotRanking <- reactive({
#      ggplot(rankingGenorm(),
#             aes(x = Rank, y = Avg.M, color = Genes)) +
#        geom_line(color = input$color) +
#        geom_point() +
#        lims(x = rankingGenorm()$Gene, y = c(0, max(rankingGenorm()$Avg.M))) +
#        ylab("Average expression stability M") +
#        xlab("<==== Most Stable Gene   Least Stable Gene ====>") +
#        ggtitle("Gene stability measure by ctrlGene Package (geNorm)") +
#        theme(
#          plot.title = element_text(size = 14, hjust = 0.5)) +
#        scale_fill_brewer(palette = "RdBu") +
#        theme_minimal()
#    })
# 
# 
#    pairwisePlotGenorm <- function(){
#      #ctvalue_new <- ctvalues()[,-length(ctvalues())]
#      pairwise <- pairwiseV(qvalue()[,2:ncol(qvalue())], ctVal = FALSE)
#      colnames(pairwise) <- c("PairWise", "Values")
#      mypalette <- brewer.pal(8, "Spectral")
#      barplot(cbind(pairwise$Values), beside = TRUE, col = mypalette, space = c(0, 0.5), ylim = c(0, 0.25), xlab = c("PairWise Variations"))
#      axis(1, at = 1:length(pairwise$Values), labels = as.character(pairwise$PairWise))
#      abline(h = seq(0, max(pairwise$Values), by = 0.2), lty = 2, col = "grey")
#      abline(h = 0.15, lty = 1, col = input$color)
#      box()
#    }
# 
# 
#    output$plotGenormShow <- renderPlot({
#      if (input$genormplot == "1") {
#        gnormPlotRanking()
#      } else if (input$genormplot == "2") {
#        pairwisePlotGenorm()
#      }
# 
#      #gnormPlotRanking()
#    })
# 
# 
#    output$downloadgenormplots <- downloadHandler(
#        filename = function() {
#          if (input$genormplot == "1") {
#            paste("genormStabilityPlot", input$type3, sep = '.')
#          } else if (input$genormplot == "2") {
#            paste("genormPairWisePlot", input$type3, sep = '.')
#          }
#        },
#        content = function(file) {
#          if (input$genormplot == "1") {
#          ggsave(file, plot = gnormPlotRanking(), device = input$type3)
#          } else if (input$genormplot == "2") {
#            if (input$type3 == "png")
#              png(file)
#             else
#              pdf(file)
#              pairwisePlotGenorm()
#              dev.off()
# 
#        }
#          })
# 
# 
#    output$downloadpairwise <- downloadHandler(
#      filename = function() {
#        paste("pairWisePlot", input$type4, sep = ".")
#      },
#      content = function(file) {
#        if (input$type4 == 'png')
#          png(file)
#        else
#          pdf(file)
#        pairwisePlotGenorm()
#        dev.off()
#      }
#    )
# 
# ################# NORMFINDER #################
# 
# ###### normfind script block
# 
# 
# 
#    normTrans <- reactive({
#      inFile1 <- input$file1
#      if (is.null(inFile1)) {return(NULL)}
#      ctvalues2 <- read.table(inFile1$datapath, header = input$header, sep = input$sep, row.names = 1)
#      ctvalueTrans <- t(ctvalues2)
# 
#        return(ctvalueTrans)
# 
#    })
# 
#    output$transtable <- renderTable({
#        normTrans()
#    },  rownames = TRUE)
# 
#    output$downloadTransNorm <- downloadHandler(
#        filename = function() {
#            paste("dataTransNormFinder", "txt", sep = '.')
#        },
#        content = function(file) {
#            write.table(normTrans(), file)
#        }
#    )
# 
#    normFinderFunction <- reactive({
#         inFile3 <- input$file3
#        if (is.null(inFile3)) {return(NULL)}
#        # ctvalues2 <- read.table(inFile3$datapath, header = input$header, sep = input$sep)
# 
# 
#        resulttotal <- Normfinder(inFile3$datapath, Groups = input$group, ctVal = TRUE)
#        return(resulttotal)
# 
# 
#    }) # script from normfinder algorithm
# 
# 
#    stabilityNofmFinder <- reactive({
#      if (input$group == TRUE) {
#        stabilityG <- data.table(rank = 1:nrow(normFinderFunction()$Ordered), normFinderFunction()$Ordered, keep.rownames = TRUE)
#        colnames(stabilityG) <- c("Rank", "Gene", "GroupDif", "GroupSD", "Stability")
#        return(stabilityG)
#      } else if (input$group == FALSE) {
#        stabilityNG <- data.table(rank = 1:nrow(normFinderFunction()$Ordered), Gene = rownames(normFinderFunction()$Ordered), normFinderFunction()$Ordered$GroupSD)
#        colnames(stabilityNG) <- c("Rank", "Gene", "Stability")
#        return(stabilityNG)
#        }
#    })
# 
#    #
#    # stabilityNofmFinderNG <- reactive({
#    #   stabilityNG <- data.table(rank = 1:nrow(normFinderFunction()$Ordered), Gene = rownames(normFinderFunction()$Ordered), normFinderFunction()$Ordered$GroupSD)
#    #   colnames(stabilityNG) <- c("Rank", "Gene", "Stability")
#    #   return(stabilityNG)
#    # })
# 
# 
#    output$normFinderTable <- renderTable({
#      if (!is.null(input$file3$datapath)) {
#        if (input$group == TRUE) {
#          if (input$normfindertable == "1") {
#             stabilityNofmFinder()
#           } else if (input$normfindertable == "2") {
#              normFinderFunction()$PairOfGenes
#         }
#           } else if (input$group == FALSE) {
#             if (input$normfindertable == "2") {
#               return(paste("This table will be generated only if there is a groups column in dataset"))
#             } else {
#               stabilityNofmFinder()
#             }
# 
#         #return(paste("Groups igual FALSE"))
#           }
#        } else {
#        return(paste("Please, download the file 'dataTransNormFinder.txt' on the previous tab and upload it again"))
#      }
#    })
# 
# 
#    output$dowRankNormFinder <- downloadHandler(
#      filename = function() {
#        paste("rankNormfinder", "csv", sep = ".")
#      },
#      content = function(file) {
#        fwrite(stabilityNofmFinder(), file, quote = FALSE, sep = "\t")
#      }
#    )
# 
#     normfinderPlot <- reactive({
#       ggplot(stabilityNofmFinder(), aes(x = Rank, y = Stability, color = Gene)) +
#         geom_line(color = input$color2) +
#         geom_point() +
#         lims(x = stabilityNofmFinder()$Gene, y = c(0, max(stabilityNofmFinder()$Stability))) +
#         ylab("Average expression stability (S)") +
#         xlab("<==== Most Stable Gene   Least Stable Gene ====>") +
#         ggtitle("Gene Stability Measure by NormFinder") +
#         theme(
#           plot.title = element_text(size = 14, hjust = 0.5)) +
#         scale_fill_brewer(palette = "RdBu") +
#         theme_minimal()
# 
#     })
# 
#    output$normFinderPlot <- renderPlot({
#      if (!is.null(input$file3$datapath)) {
#        normfinderPlot()
#      } else {
#        return(paste("Please, download the file 'dataTransNormFinder.txt' on the previous tab and upload it again"))
#      }
#    })
# 
#    output$downnormfinderplot <- downloadHandler(
#      filename = function() {
#        paste("normFinderStabilityPlot", input$type5, sep = ".")
#      },
#      content = function(file) {
#        ggsave(file, normfinderPlot(), device = input$type5)
#      }
#    )
# 
# 
# 
# ########################### BESTKEEPER ############################
# 
#    bestkeepanaly <- reactive({
#      inFile1 <- input$file1
#      if (is.null(inFile1)) {return(NULL)}
#      ctvalue <- read.table(inFile1$datapath, header = input$header, sep = input$sep, row.names = 1)
#      if (input$group == TRUE) {
#      bestkeeper_results = bestKeeper(ctvalue[,1:ncol(ctvalue) - 1], ctVal = TRUE)
#      return(bestkeeper_results)
#      } else if (input$group == FALSE) {
#        bestkeeper_results = bestKeeper(ctvalue[,1:ncol(ctvalue)], ctVal = TRUE)
#        return(bestkeeper_results)
#      }
#    })
# 
#    cpstatbestkeeper <- reactive({
#      cpStatBestkeeper <- data.table(bestkeepanaly()$CP.statistics, keep.rownames = TRUE)
#      return(cpStatBestkeeper)
#    })
# 
#    pair.Wise.corbestkeeper <- reactive({
#      pairwisebestkeeper <- data.table(bestkeepanaly()$pair.Wise.cor, keep.rownames = TRUE)
#      return(pairwisebestkeeper)
#    })
# 
#    HKG.vs.BestKeeper <- reactive({
#      hkgvsbestkeeper <- data.table(bestkeepanaly()$HKG.vs.BestKeeper, keep.rownames = TRUE)
#      return(hkgvsbestkeeper)
#    })
# 
#    finalbestrank <- reactive({
#      bestkeeper_genes <- t(sort(bestkeepanaly()$CP.statistics[6,]));
#      best_genes <- t(bestkeeper_genes);
#      best_genes_values <- data.table(Rank = 1:nrow(best_genes), rownames(best_genes), best_genes)
#      colnames(best_genes_values) <- c('Rank', 'Genes', 'SD_Value')
#      return(best_genes_values)
#    })
# 
#     output$bestkeepertables <- renderTable({
#       if (input$table == "1") {
#         cpstatbestkeeper()
#         #as.data.frame(bestkeepanaly()$CP.statistics)
#         }
#       else if (input$table == "2") {
#         pair.Wise.corbestkeeper()
#         #bestkeepanaly()$pair.Wise.cor
#         }
#       else if (input$table == "3") {
#         HKG.vs.BestKeeper()
#         #bestkeepanaly()$HKG.vs.BestKeeper
#       }
#       else if (input$table == "4") {
#         finalbestrank()
# }
# })
# 
# 
#     output$dowbestkeeper <- downloadHandler(
#       filename = function() {
#         if (input$table == "1") {
#           paste("CPstatisticsBestkeeper", "csv", sep = ".")
#         }
#         else if (input$table == "2") {
#           paste("pairWisecorBestkeeper", "csv", sep = ".")
#         }
#         else if (input$table == "3") {
#           paste("HKG.vs.BestKeeperNestkeeper", "csv", sep = ".")
#         }
#         else if (input$table == "4") {
#         paste("finalRankingBestkeeper", "csv", sep = ".")
#         }
#       }
#         ,
#       content = function(file) {
#         if (input$table == "1") {
#           fwrite(cpstatbestkeeper(), file, quote = FALSE)
#         }
#         else if (input$table == "2") {
#           fwrite(pair.Wise.corbestkeeper(), file, quote = FALSE, sep = "\t")
#         }
#         else if (input$table == "3") {
#           fwrite(HKG.vs.BestKeeper(), file, quote = FALSE, sep = "\t")
#         }
#         else if (input$table == "4") {
#           fwrite(finalbestrank(), file, quote = FALSE, sep = "\t")
#         }
#       }
#     )
# 
# 
#     plotreacbestkeepr <- reactive({
#       ggplot(finalbestrank(),
#              aes(x = Rank, y = SD_Value, color = Genes)) +
#         geom_line(color = input$color3) +
#         geom_point() +
#         lims(x = finalbestrank()$Genes, y = c(0, max(finalbestrank()$SD_Value))) +
#         ylab("Average expression stability (Standard Deviation - SD)") +
#         xlab("<==== Most Stable Gene   Least Stable Gene ====>") +
#         ggtitle("Gene Stability Measure by Bestkeeper") +
#         theme(
#           plot.title = element_text(size = 14, hjust = 0.5)) +
#         scale_fill_brewer(palette = "RdBu") +
#         theme_minimal()
#     })
# 
# 
#     output$plotbestkeeper <- renderPlot({
#       plotreacbestkeepr()
#     })
# 
#     output$dowFinalRankBestkeeper <- downloadHandler(
#       filename = function() {
#         paste("bestkeeperPlot", input$type6, sep = ".")
#       },
#       content = function(file) {
#         ggsave(file, plotreacbestkeepr(), device = input$type6)
#       }
#     )
# 
# 
# 
# 
#     CESP <- reactive({
# 
#       genes <- strsplit(rankingGenorm()$Genes[1], split = "-")
#       # gene1 <- substr(basename(paste0(rankingGenorm()$Genes[1])), 1, nchar(basename(paste0(rankingGenorm()$Genes[1]))) - 7)
#       # gene2 <- substr(basename(paste0(rankingGenorm()$Genes[1])), 8, nchar(basename(paste0(rankingGenorm()$Genes[1]))) - 0)
# 
#       ratools <- t(data.frame(NormFinder = stabilityNofmFinder()$Gene, GeNorm = c(genes[[1]][1], genes[[1]][2], rankingGenorm()$Genes[2:length(rankingGenorm()$Genes)]), Bestkeeper = finalbestrank()$Genes))
#       colnames(ratools) <- paste(1:ncol(ratools))
# 
# 
#       ra_weight_ratools <- t(data.frame(NormFinder = stabilityNofmFinder()$Stability, GeNorm = c(rankingGenorm()$Avg.M[1], rankingGenorm()$Avg.M), Bestkeeper = finalbestrank()$SD_Value))
#       colnames(ra_weight_ratools) <- paste(1:ncol(ratools))
# 
# 
#       k <- ncol(ratools)
# 
#       #method = opt$method
#       #distance = opt$distance
#       #maxIter = opt$iteraction
# 
# 
#       if (k >= 10) {
#         CESP <- RankAggreg(ratools, k = ncol(ratools), ra_weight_ratools, method = input$method, distance = input$distance, weight = .25, rho = .1, maxIter = 1000, verbose = TRUE)
#       } else {
#         CESP <- BruteAggreg(ratools, k = ncol(ratools), ra_weight_ratools, distance = input$distance)
#         }
#       return(CESP)
#     })
# 
# 
# 
#     rankfinal <- reactive({
# 
#       final_res <- data.frame(Rank = paste(1:length(CESP()$top.list)), Gene = CESP()$top.list)
#         return(final_res)
# 
#     })
# 
# 
#     output$finalranking <- renderTable({
#       rankfinal()
# 
#     })
# 
# 
#     output$dowfinalrank <- downloadHandler(
#       filename = function() {
#         paste("finalRanking", "csv", sep = ".")
#       },
#       content = function(file) {
#         fwrite(rankfinal(), file, quote = FALSE, sep = "\t")
#       }
#     )
# 
# 
#       output$plotfinal <- renderPlot({
# 
#         plot(CESP())
#         })
# 
# 
#     output$dowfinalPlot <- downloadHandler(
#       filename = function() {
#         paste("finalRankingPlot", input$type7, sep = ".")
#       },
#       content = function(file) {
#         if (input$type7 == 'png')
#           png(file)
#         else
#           pdf(file)
#           plot(CESP())
#           dev.off()
# 
#       }
#     )
# 
# ######### relative expression analysis
# 
#     targetGenes <- reactive({
#       inFile4 <- input$file4
#       if (is.null(inFile4)) {return(NULL)}
#       trget <- fread(inFile4$datapath, header = input$header, sep = input$sep, data.table = FALSE)
# 
#       #colnames(dataset) <- c("SampleName", colnames(dataset[,2:ncol(dataset)]))
#       return(trget)
#     })
# 
# 
#     efficiencyTarget <- reactive({
#       inFile5 <- input$file5
#       if (is.null(inFile5)) {return(NULL)}
#       effTarg <- fread(inFile5$datapath, header = input$header2, sep = input$sep2, col.names = c('GENE', 'EFFICIENCY'))
# 
#       return(effTarg)
#     })
# 
# 
#    geom_mean <- reactive({
#       gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
#         if (any(x < 0, na.rm = TRUE)) {
#           return(NaN)
#         }
#         if (zero.propagate) {
#           if (any(x == 0, na.rm = TRUE)) {
#             return(0)
#           }
#           exp(mean(log(x), na.rm = na.rm))
#         } else {
#           exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
#         }
#       }
#       #input$normQuant
#       ctvalue_new <- ctvalues()[,-length(ctvalues())]
#       bestGenes <- select(.data = ctvalue_new, rankfinal()[1:input$normQuant,2])
#       geomean <- array()
#       for (i in 1:nrow(bestGenes)) {
#         geomean[i] <- gm_mean(bestGenes[i,])
#       }
#       return(geomean)
#     })
# 
#    bestGenes <- reactive({
#      ctvalue_new <- ctvalues()[,-length(ctvalues())]
#      bestGenes <- select(.data = ctvalue_new, rankfinal()[1:input$normQuant,2])
#      return(bestGenes)
#    })
# 
# 
#    finalDataaa <- reactive({
#      if (is.null(input$file4$datapath)) {
#        finaldat <- data.frame(SampleName = ctvalues()$SampleName, bestGenes(), "Groups" = ctvalues()[,ncol(ctvalues())])
#        return(finaldat)
#      } else {
#        finaldat <- data.frame(SampleName = ctvalues()$SampleName, targetGenes(), bestGenes(), "Groups" = ctvalues()[,ncol(ctvalues())])
#        return(finaldat)
#      }
#    })
# 
# 
# 
#     showGenes <- reactive({
#       bestGenes <- rankfinal()[1:input$normQuant,2]
#       return(bestGenes)
#     })
# 
#     output$usedGenes <- renderText({
#       if (!is.null(input$file4$datapath) & input$data == 4) {
#         if (input$Normethod == "Vandesompele") {
#           return(c("Using", showGenes(), " as Reference Genes"))
#         }
#         else {
#           return(c("Using", showGenes(), "genes as Normalizer Factor by geometric mean"))
#         }
#       }
#     })
# 
#     normalizationMethod <- reactive({
#       bestGenes <- select(.data = ctvalues(), rankfinal()[1:input$normQuant,2])
#       datagrid <- data.frame(targetGenes(), "NormFactor" = geom_mean())
#       normFactor <- data.table("NormFactor" = geom_mean())
#       if (input$group) {
#         Groups <- ctvalues()[,ncol(ctvalues())]
#       }
#       if (input$Normethod == "Vandesompele") {
#         ctvalue_new <- ctvalues()[,-length(ctvalues())]
#         datagrid <- data.frame(targetGenes(), bestGenes)
#         n <- ncol(targetGenes()) + 1
#         if (ncol(targetGenes()) > 1) {
#           normdata <- normPCR(datagrid, c(n:ncol(datagrid)))
#           if (input$group) {
#             finaldatagrid <- data.frame(SampleName = ctvalues()$SampleName, normdata, "Groups" = ctvalues()[,ncol(ctvalues())])
#             return(finaldatagrid)
#           } else {
#             finaldatagrid <- data.frame(SampleName = ctvalues()$SampleName, normdata)
#             return(finaldatagrid)
#           }
#         } else {
#           normdata <- normPCR(datagrid, c(n:ncol(datagrid)))
#           normdatatable <- data.table(normdata)
#           colnames(normdatatable) <- names(targetGenes())
#           if (input$group) {
#             finaldatagrid <- data.frame(SampleName = ctvalues()$SampleName, normdatatable, "Groups" = ctvalues()[,ncol(ctvalues())])
#             return(finaldatagrid)
#           } else {
#             finaldatagrid <- data.frame(SampleName = ctvalues()$SampleName, normdatatable)
#             return(finaldatagrid)
#           }
#         }
#       }
#       else if (input$Normethod == "Livak") {
#         if (ncol(targetGenes()) > 1) {
#           dctData <- mutate_at(.tbl = datagrid, vars(-matches("NormFactor")), list(dCt = ~ . - NormFactor)) %>%
#             rename_at(vars(contains("_dCt") ), list(~paste("(∆Ct)", gsub("_dCt", "", .), sep = "")))
#         } else {
#           dctData <- targetGenes()[,] - normFactor[,]
#           colnames(dctData) <- paste0("(∆Ct)", names(targetGenes()))
#         }
# 
#         if (input$group) {
#             dctDataControl <- data.table(dctData, "Groups" = Groups) %>%
#               filter(Groups == input$Contaverage) %>%
#               select(starts_with("(∆Ct)"))
# 
#             controlAverage <- colMeans(dctDataControl)
# 
#             dctDatatargets <- data.table(dctData, "Groups" = Groups) %>%
#               select(starts_with("(∆Ct)"))
# 
#             listContAver <- list()
#             for (i in 1:length(controlAverage)) {
#               listContAver[i] <- controlAverage[i]
#             }
# 
#         } else {
#             dctDataControl <- data.table(dctData) %>%
#               select(starts_with("(∆Ct)"))
# 
#             controlAverage <- apply(dctDataControl, 2, min)
# 
#             dctDatatargets <- data.table(dctData) %>%
#               select(starts_with("(∆Ct)"))
# 
#             listContAver <- list()
#             for (i in 1:length(controlAverage)) {
#               listContAver[i] <- controlAverage[i]
#             }
#             listContAver <- do.call('cbind', listContAver)
# 
#         }
#         ddct <- dctDatatargets - listContAver
#         ddct <- rename_at(ddct, vars(contains("(∆Ct)")), list(~paste(gsub("(∆Ct)", "∆∆Ct", .), sep = "")))
# 
# 
# 
#           if (!is.null(efficiencyTarget())) {
#             # nE <- nrow(efficiencyTarget()) - input$normQuant
#             # neweff <- efficiencyTarget()[1:nE,]
#             listeff <- list()
#             for (i in 1:nrow(efficiencyTarget()[,2])) {
#               listeff[i] <- efficiencyTarget()[i,2]
#             }
#             doisddct <- listeff^-(ddct[,1:ncol(ddct)])
#             doisddct <- rename_at(doisddct, vars(contains("(∆∆Ct)")), list(~paste("2^-", gsub("∆∆Ct_", "()", .), sep = "")))
#           } else {
#             doisddct <- 2^-(ddct[,])
#             doisddct <- rename_at(doisddct, vars(contains("(∆∆Ct)")), list(~paste("2^-", gsub("∆∆Ct_", "()", .), sep = "")))
#           }
# 
# 
#           if (input$group) {
#             finaldatagrid <- data.table(SampleName = ctvalues()$SampleName, dctData, ddct, doisddct, Groups)
#             return(finaldatagrid)
#           } else {
#             finaldatagrid <- data.table(SampleName = ctvalues()$SampleName, dctData, ddct, doisddct)
#             return(finaldatagrid)
#           }
#       }
#       else if (input$Normethod == "Pfaffl") {
#         if (!is.null(targetGenes()) & !is.null(efficiencyTarget())) {
#           if (nrow(efficiencyTarget()) == ncol(targetGenes())) {
# 
#             inFile2 <- input$file2
#             inFile4 <- input$file4
#             inFile5 <- input$file5
# 
#             #refereceGenes <- bestGenes
# 
#             efficiencyRG <- fread(inFile2$datapath, header = input$header2, sep = input$sep2, col.names = c('GENE', 'EFFICIENCY'), data.table = FALSE)
# 
#             target <- fread(inFile4$datapath, header = input$header, sep = input$sep, data.table = FALSE)
# 
#             efficiencyTG <- fread(inFile5$datapath, header = input$header2, sep = input$sep2, col.names = c('GENE', 'EFFICIENCY'))
# 
# 
# 
#             normFactor <- data.table("NormFactor" = geom_mean())
#             if (input$group) {
#               Groups <- ctvalues()[,ncol(ctvalues())]
#               fulldata <- data.frame(target, normFactor, "Groups" = Groups) %>%
#                 filter(Groups == input$Contaverage)
# 
#               fulldata <- fulldata[,-length(fulldata)]
# 
#               controlAverage <- colMeans(fulldata)
# 
#               dctData <- data.table(targetGenes(), normFactor)
# 
#               listContAver <- list()
#               for (i in 1:length(controlAverage)) {
#                 listContAver[i] <- controlAverage[i]
#               }
#             } else {
#               fulldata <- data.frame(target, normFactor)
#               #fulldata <- fulldata[,-length(fulldata)]
#               controlAverage <- apply(fulldata, 2, min)
# 
#               dctData <- data.table(targetGenes(), normFactor)
# 
#               listContAver <- list()
#               for (i in 1:length(controlAverage)) {
#                 listContAver[i] <- controlAverage[i]
#               }
#             }
# 
#             dct <- data.frame(dctData - listContAver)
#             #dct <- data.frame(dct1)
#             colnames(dct) <- paste0("dCt_", names(dctData))
# 
# 
# 
#             if (ncol(target) == 1) {
#               dctTarget <- data.table(dct[,-length(dct)])
#               dctRG <- data.table(dct[,ncol(dct)])
#               colnames(dctTarget) <- names(target)
#               colnames(dctRG) <- names(normFactor)
#             } else {
#               dctTarget <- data.table(dct[,-length(dct)])
#               dctRG <- data.table(dct[,ncol(dct)])
#             }
# 
#             listEffRG <- list()
#             for (i in colnames(bestGenes)) {
#               listEffRG[[i]] <- efficiencyRG[efficiencyRG[, 1] == i,]
#             }
# 
#             effRG <- do.call("rbind", listEffRG)
# 
#             #effRG <- as.data.frame()
#             effRG <- mean(effRG$EFFICIENCY)
# 
# 
#             # # nE <- nrow(efficiency()) - normQuant
#             # # neweff <- efficiency()[1:nE,]
#             #
#             #
#             effTG <- list()
#             for (i in 1:nrow(efficiencyTG)) {
#               effTG[i] <- efficiencyTG[i,2]
#             }
# 
# 
# 
#             EdctTarget <- data.frame(effTG^dctTarget[,])
#             EdctRG <- effRG^dctRG
# 
#             expratio <- list()
#             for (i in 1:ncol(EdctTarget)) {
#               expratio[[i]] <- EdctTarget[,i]/EdctRG
#             }
# 
#             expratioData <- data.frame(expratio)
#             colnames(expratioData) <- paste0("ExpressionRatio_", names(target))
#             if (input$group) {
#               finaldatagrid <- data.frame(SampleName = ctvalues()$SampleName, dct, expratioData, Groups)
#             } else {
#             finaldatagrid <- data.frame(SampleName = ctvalues()$SampleName, dct, expratioData)
#             }
#             return(finaldatagrid)
# 
#           } else {
#             print("You need to give the efficiency value for the Reference Gene (s) for running the pffafl method.")
#           }
#         } else {
#           return("Please upload the Ct and the efficiency files")
#         }
#       }
#     })
# 
# 
#     output$normTables <- renderTable({
#       normalizationMethod()
#     })
# 
#     output$legendTables <- renderText({
#       if (!is.null(input$file4$datapath) & input$data == 4) {
#         if (input$Normethod == "Vandesompele") {
#           #return(c("NormFactor: Geometric mean of the most stable reference genes"))
#         }
#         else if (input$Normethod == "Livak") {
#           gsub(pattern = "\\n", replacement = "<br/>", paste("NormFactor: Geometric mean of the most stable reference genes;", "∆Ct: Delta Ct of the Target Gene;", "∆∆Ct: DeltaDelta Ct of the Target Gene;", "2^-(∆∆Ct): Fold gene expression values;", sep = "\n"))
#           #paste("NormFactor: Geometric mean of the most stable reference genes;", "∆Ct: Delta Ct of the Target Gene;", "∆∆Ct: DeltaDelta Ct of the Target Gene;", "2^-(∆∆Ct): Fold gene expression values;", sep = "\n")
#         }
#         else if (input$Normethod == "Pfaffl") {
#           #return(c("Using", showGenes(), "genes as Normalizer Factor by geometric mean"))
#         }
#       }
#     })
#   
#     
#     
#     output$targetTable <- renderTable({
#       if (is.null(input$file4$datapath)) {
#         if (input$headss == TRUE) {
#           finalDataaa()
#         } else {
#           head(finalDataaa())
#         }
#       } else {
#         if (input$data == 1) {
#           if (input$headss == TRUE) {
#             data.table(SampleName = ctvalues()$SampleName, targetGenes(), "Groups" = ctvalues()[,ncol(ctvalues())])
#           } else {
#             head(data.table(SampleName = ctvalues()$SampleName, targetGenes(), "Groups" = ctvalues()[,ncol(ctvalues())]))
#           }
#         }
#         else if (input$data == 2) {
#           if (input$headss == TRUE) {
#             efficiencyTarget()
#           } else {
#             head(efficiencyTarget())
#           }
#         }
#         else if (input$data == 3) {
#           if (input$headss == TRUE) {
#             finalDataaa()
#           } else {
#             head(finalDataaa())
#           }
#         }
#         else if (input$data == 4) {
#           if (!is.null(input$file4$datapath)) {
#             if (input$headss == TRUE) {
#               normalizationMethod()
#             } else {
#               head(normalizationMethod())
#             }
#           } else {
#             return("Please upload the target genes dataset")
#           }
#         }
#       }
#     })
#     
# 
#     output$downRelExp <- downloadHandler(
#       filename = function() {
#         if (input$Normethod == "Vandesompele") {
#           paste("dataNormalizedVandesompele", "csv", sep = ".")
#         }
#         else if (input$Normethod == "Livak") {
#           paste("dataNormalizedLivak", "csv", sep = ".")
#         }
#         else if (input$Normethod == "Pfaffl") {
#           paste("dataNormalizedPfaffl", "csv", sep = ".")
#         }
#       }, 
#       content = function(file) {
#         fwrite(normalizationMethod(), file, quote = FALSE, sep = "\t")
#       }
#     )
#     
#     
#     
#     
#     
#     vandesompelenExplanation <- a("Click here", href = "https://toptipbio.com/qpcr-multiple-reference-genes/", target = "_blank")
#     livakExplanation <- a("Click here", href = "https://toptipbio.com/delta-delta-ct-pcr/", target = "_blank")
#     pfaffilExplanation <- a("Click here", href = "https://toptipbio.com/pfaffl-method-qpcr/#:~:text=The%20Pfaffl%20method%2C%20named%20after,Nucleic%20Acids%20Research%20in%202001.", target = "_blank")
#     
#     url_RGct <- a("Click here", href = "https://github.com/hanielcedraz/enDogenesShiny/blob/main/endogenesCtFile.txt", target = "_blank")
#     
#     url_eff <- a("Click here", href = "https://github.com/hanielcedraz/enDogenesShiny/blob/main/endogenesEfficiencyFile.txt", target = "_blank")
#     
#     output$example_ct <- renderUI({
#       tagList(url_RGct, "to download an example of the txt file with ct data.")
#     })
#     
#     output$example_eff <- renderUI({
#       tagList(url_eff, "to download an example of the txt file with efficiency data.")
#     })
# 
#     RankAggreg <- a("RankAggreg, ", href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-62", target = "_blank")
#     ctrlGene <- a("ctrlGene, ", href = "https://cran.r-project.org/web/packages/ctrlGene/index.html", target = "_blank")
#     SLqPCR <- a("SLqPCR, ", href = "https://www.bioconductor.org/packages/release/bioc/html/SLqPCR.html")
#     normfinder <- a("normfinder", href = "https://moma.dk/normfinder-software", target = "_blank")
#     
#     output$links <- renderUI({
#       tagList("This shinyApp uses the following packages: ", RankAggreg, ctrlGene, SLqPCR, "and, a function for", normfinder, "provided by MOMA -Department of Molecular Medicine")
#     })
#     
#     
#     output$analysisExplanation <- renderUI({
#       if (input$Normethod == "Vandesompele") {
#         tagList(vandesompelenExplanation, "for more details about the Vandesompele Method")
#       }
#       else if (input$Normethod == "Livak") {
#         tagList(livakExplanation, "for more details about the Livak Method")
#       }
#       else if (input$Normethod == "Pfaffl") {
#         tagList(pfaffilExplanation, "for more details about the Pfaffl Method")
#       }
#     })
#     
#     url_TGct <- a("Click here", href = "https://github.com/hanielcedraz/enDogenesShiny/blob/main/targetsCtFile.txt", target = "_blank")
#     
#     output$targetCtExample <- renderUI({
#       tagList(url_TGct, "to download an example of the txt file with ct data.")
#     })
    
    output$tb <- renderUI({
      fluidRow(
        tabBox(
          id = "tabs",
          title = "",
          selected = "PendingOrders",
          status = "primary",
          solidHeader = FALSE, 
          type = "tabs",
          width = 12,
          tabPanel(
            title = "DataSets",
            mainPanel(
              width = 12,
              selectInput(
                "dataset", "Select the Table",
                choices = c("Endogenous Dataset" = "1",
                            "Normalized Data" = "2",
                            "Efficiency File" = "3"
                )
              ),
              tableOutput(
                outputId = "tables"
              ),
              tipify(
                checkboxInput(
                  inputId = 'head', 
                  label = "Expand", 
                  value = FALSE
                ),
                "Click here to show the entire dataset",
                placement = "left",
                trigger = "hover"
              ),
              downloadButton(
                outputId = "downloadNormData", 
                label = "Download the Norm data"
              )
            ),
          ),
          tabPanel(
            title = "Boxplot ct Values",
            mainPanel(
              width = 12,
              box(
                title = "MainTable",
                id = "box1",
                selected = "",
                status = "primary",
                solidHeader = FALSE,
                #type = "tabs",
                width = 12,
                side = "right",
                selectInput(
                  inputId = "plot", "Select the Plot",
                  choices = c("Boxplot by Gene" = "gene",
                              "Boxplot by Group" = "group"
                  )
                ),
                plotOutput(
                  outputId = "boxplots"
                ),
                downloadButton(
                  outputId = "downloadBoxPlots",
                  label = "Download the Plot"
                ),
                radioButtons(
                  inputId = "type", 
                  label = "Choose the file type",
                  choices = list('png', 'pdf', "tiff"),
                  selected = 'png'
                )
              )
            ),
          ),
          tabPanel(
            title = "Ct values transposed",
            mainPanel(
              width = 12,
              h5(helpText("This dataset will be used
                         in the next analysis. Please, 
                           download it")),
              tableOutput(
                outputId = "transtable"
              ),
              downloadButton(
                outputId = "downloadTransNorm", 
                label = "Download the transposed dataset"
              ),
              fileInput(inputId = "file3", 
                        label = "Please, download the above dataset 
                              and upload it again in this box. This dataset 
                              will be used in the next analysis"
              )
            )
          ),
          tabPanel(
            title = "NormFinder Results",
            mainPanel(
              width = 12,
              br(),
              uiOutput(
                outputId = "downloadHistory"
              ),
              br(),
              dataTableOutput(
                outputId = "tableHistory"
              )
            )
          ),
          tabPanel(
            title = "NotifyFiles",
            mainPanel(
              width = 12,
              box(
                selectInput(
                  inputId = "normfindertable", 
                  label = "Select the Table",
                  choices = c("NormFinder Ranking" = "1",
                              "Pair of Genes" = "2"
                  )
                ),
                tableOutput(
                  outputId = "normFinderTable"
                ),
                downloadButton(
                  outputId = "dowRankNormFinder", 
                  label = "Download the file"
                )
              ),
              box(
                plotOutput(
                  outputId = "normFinderPlot"
                ),
                downloadButton(
                  outputId = "downnormfinderplot", 
                  label = "Download the file"
                ),
                textInput(
                  inputId = "color2", 
                  label = "Choose the line color",
                  value = "black", 
                  placeholder = "Type a valid color"
                ),
                radioButtons(
                  inputId = "type5", 
                  label = "Choose the file type",
                  choices = list('png', 'pdf', "tiff"),
                  selected = 'png'
                )
              )
            )
          ),
          tabPanel(
            title = "geNorm Results",
            mainPanel(
              width = 12,
              box(
                selectInput(
                  inputId = "genormtable", 
                  label = "Select the Table",
                  choices = c("Genorm Ranking" = "1",
                              "PairWise Correlation" = "2"
                  )
                ),
                tableOutput(
                  outputId = "geNormtable"
                ),
                downloadButton(
                  outputId = "downloadGeNormranking",
                  label = "Download the geNorm Ranking"
                )
              ),
              box(
                selectInput(
                  inputId = "genormplot", 
                  label = "Select the Plot",
                  choices = c("Genorm Ranking" = "1",
                              "PairWise Correlation" = "2"
                  )
                ),
                plotOutput(
                  outputId = "plotGenormShow"
                ),
                downloadButton(
                  outputId = "downloadgenormplots",
                  label = "Download the Plot"
                ),
                textInput(
                  inputId = "color", 
                  label = "Choose the line color", 
                  value = "black", 
                  placeholder = "Type a valid color"
                ),
                radioButtons(
                  inputId = "type3", 
                  label = "Choose the file type",
                  choices = list('png', 'pdf', "tiff"),
                  selected = 'png'
                )
              )
            )
          ),
          tabPanel(
            title = "Bestkeeper Results",
            mainPanel(
              width = 12,
              selectInput(
                inputId = "table", 
                label = "Select the table",
                choices = c("Statistic Summary" = "1",
                            "PairWise Correlation" = "2",
                            "HKG vs. Bestkeeper" = "3",
                            "Bestkeeper Ranking" = "4"
                )
              ),
              tableOutput(
                outputId = "bestkeepertables"
              ),
              downloadButton(
                outputId = "dowbestkeeper", 
                label = "Download the file"
              ),
              plotOutput(
                outputId = "plotbestkeeper"
              ),
              textInput(
                inputId = "color3", 
                label = "Choose the line color", 
                value = "black", 
                placeholder = "Type a valid color"
              ),
              downloadButton(
                outputId = "dowFinalRankBestkeeper", 
                label = "Download the file"
              ),
              radioButtons(
                inputId = "type6", 
                label = "Choose the file type",
                choices = list('png', 'pdf', "tiff"),
                selected = 'png'
              )
            )
          ),
          tabPanel(
            title = "Final Ranking",
            mainPanel(
              width = 12,
              sh5(helpText("Please wait until the results to be generate. It could take some time")),
              tableOutput(
                outputId = "finalranking"
              ),
              tipify(
                radioButtons(
                  inputId = "method", 
                  label = "Choose the method",
                  choices = list("CE", "GA"),
                  selected = "CE", 
                  inline = TRUE
                ),
                title = "Cross Entropy Monte Carlo (CE) or Genetic Algorithm (GA)",
                placement = "left", 
                trigger = "hover"
              ),
              tipify(
                radioButtons(
                  inputId = "distance", 
                  label = "Choose the distance",
                  choices = list("Spearman", "Kendall"),
                  selected = "Spearman", 
                  inline = TRUE
                ),
                title = "Distance to be used in the similarity",
                placement = "left", 
                trigger = "hover"
              ),
              downloadButton(
                outputId = "dowfinalrank", 
                label = "Download the file"
              ),
              plotOutput(
                outputId = "plotfinal"
              ),
              #h5(helpText("Please wait until the plot to be generate. It could take some time")),
              downloadButton(
                outputId = "dowfinalPlot", 
                label = "Download the file"
              ),
              radioButtons(
                inputId = "type7", 
                label = "Choose the file type",
                choices = list('png', 'pdf', "tiff"),
                selected = 'png'
              )
            )
          ),
          tabPanel(
            title = "Relative Expression",
            mainPanel(
              uiOutput(
                outputId = "targetCtExample"
              ),
              div(style = "display:inline.block",
                  splitLayout(
                    fileInput(
                      inputId = "file4", 
                      label = "Upload the Target Genes ct file"
                    ),
                    fileInput(
                      inputId = "file5", 
                      label = "Upload the Target Genes efficiency file"
                    )
                  )
              ),
              selectInput(
                inputId = "data", 
                label = "Select the table",
                choices = c("Target Genes Ct" = "1",
                            "Target Genes Efficiency" = "2",
                            "Target and Reference Gene(s)" = "3",
                            "Normalization Tables" = "4")
              ),
              tipify(
                radioButtons(
                  inputId = "Normethod", 
                  label = "Choose the method",
                  choices = list("Vandesompele", "Livak", "Pfaffl"),
                  selected = "Vandesompele", 
                  inline = TRUE
                ),
                title = "Choose the method used for normalization",
                placement = "left", 
                trigger = "hover"
              ),
              div(
                style = "display:inline.block",
                splitLayout(
                  textInput(
                    inputId = "normQuant", 
                    label = HTML(
                      paste("Choose the number of genes to use as reference.", "The most stable genes will be used", sep = "<br/>")
                    ), 
                    value = 3, 
                    placeholder = "Default = 3"
                  ),
                  textInput(
                    inputId = "Contaverage", 
                    label = "Choose the group to be used as Control average.", 
                    value = "Control", 
                    placeholder = "Default = Control"
                  )
                )
              ),
              textOutput(
                outputId = "usedGenes"
              ),
              tableOutput(
                outputId = "targetTable"
              ),
              htmlOutput(
                outputId = "legendTables"
              ),
              uiOutput(
                outputId = "analysisExplanation"
              ),
              tipify(
                checkboxInput(
                  inputId = 'headss', 
                  label = "Expand", 
                  value = FALSE
                ),
                title = "Click here to show the entire dataset",
                placement = "left", 
                trigger = "hover"
                ),
              downloadButton(
                outputId = "downRelExp", 
                label = "Download the table"
                )
              )
            )
          )
        )
    })
  }
)   
