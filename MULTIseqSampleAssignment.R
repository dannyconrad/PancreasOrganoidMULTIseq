###################################
# MULTI-seq Sample Assignment     #
# Pancreatic Organoids            #
# Danny Conrad                    #
# 2/15/21                         #
###################################

# Set working directory
wd <- "/Volumes/DannyData/Ling_MULTIseq_1"
setwd(wd)

# Load package
library(deMULTIplex) ## devtools::install_github('chris-mcginnis-ucsf/MULTI-seq')
source("R_scripts/classifyCells_DC.R") ## had to create a modified version of classifyCells function to enable manual thresholding for Bar11

# Make lane-specific lists of cells first (without )
cells.L1 <- colnames(L1.filt.data)
cells.L2 <- colnames(L2.filt.data)

## Define vectors for reference barcode sequences and cell IDs
bar.ref <- c("AGATTCCG", "GCAGTATC", "ATACGTGC", "CAAGGTCT", "GCATGTAC", "TTACGGTG", "AAAGGGGA", "CCATGGCG", "GGAAGGAA", "TTAGCCAG", "CGACCAGC") #MULTIseq BCs C12-C4, A11-A10
# names(bar.ref) <- c(paste("C",12:4,sep = ""), "A11", "A10")
names(bar.ref) <- c("Progenitor_1",
                    "Progenitor_2",
                    "DuctOrg8_1",
                    "DuctOrg8_3",
                    "DuctOrg16_1",
                    "DuctOrg16_2",
                    "AcinOrg8_1",
                    "AcinOrg8_2",
                    "AcinOrg8_3",
                    "AcinOrg16_1",
                    "AcinOrg16_2")

## Pre-process MULTI-seq sample barcode FASTQs
readTable_1 <- MULTIseq.preProcess(R1 = 'multi_barcodes/L1_MULTI_S3_L001_R1_001.fastq.gz', 
                                   R2 = 'multi_barcodes/L1_MULTI_S3_L001_R2_001.fastq.gz', 
                                   cellIDs = cells.L1, 
                                   cell=c(1,16), 
                                   umi=c(17,28), 
                                   tag=c(1,8))

saveRDS(readTable_1, "multi_barcodes/read.table_1.rds")

readTable_2 <- MULTIseq.preProcess(R1 = 'multi_barcodes/L2_MULTI_S4_L001_R1_001.fastq.gz', 
                                   R2 = 'multi_barcodes/L2_MULTI_S4_L001_R2_001.fastq.gz', 
                                   cellIDs = cells.L2, 
                                   cell=c(1,16), 
                                   umi=c(17,28), 
                                   tag=c(1,8))

saveRDS(readTable_2, "multi_barcodes/read.table_2.rds")


## Perform MULTI-seq sample barcode alignment
barTable_1 <- MULTIseq.align(readTable_1, cells.L1, bar.ref)

saveRDS(barTable_1, "multi_barcodes/bar.table_1.rds")

barTable_2 <- MULTIseq.align(readTable_2, cells.L2, bar.ref)

saveRDS(barTable_2, "multi_barcodes/bar.table_2.rds")


## Visualize barcode space
barTSNE_1 <- barTSNE(barTable_1[,1:11]) 
barTSNE_2 <- barTSNE(barTable_2[,1:11]) 

pdf("multi_barcodes/bc.check_1.pdf")
for (i in 3:ncol(barTSNE_1)) {
  g <- ggplot(barTSNE_1, aes(x = TSNE1, y = TSNE2, color = barTSNE_1[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(barTSNE_1)[i]) +
    theme(legend.position = "none") 
  print(g)
}
dev.off()

pdf("multi_barcodes/bc.check_2.pdf")
for (i in 3:ncol(barTSNE_2)) {
  g <- ggplot(barTSNE_2, aes(x = TSNE1, y = TSNE2, color = barTSNE_2[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(barTSNE_2)[i]) +
    theme(legend.position = "none") 
  print(g)
}
dev.off()

## Lane 1 Sample Classification------------------------------------------------------------------------------

## Round 1 
## Perform Quantile Sweep
bar.table <- barTable_1

bar.table.full <- bar.table[,1:11]
good.bars <- paste("Bar",1:11,sep="")  # NOTE: In this hypothetical example, barcodes 91-96 were not detected
bar.table <- bar.table.full[, good.bars]  # Remove missing bars and summary columns
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))

## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells_DC(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema)) # Used 1.2 as a threshold for Bar11
# round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

## Round 2
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells_DC(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema)) # Used 1.2 as a threshold for Bar11
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])

## Repeat until all no negative cells remain (usually 3 rounds)...
final.calls_1 <- c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls_1) <- c(names(round2.calls),neg.cells)


## Lane 2 Sample Classification------------------------------------------------------------------------------

## Round 1 
## Perform Quantile Sweep
bar.table <- barTable_2

bar.table.full <- bar.table[,1:11]
good.bars <- paste("Bar",1:11,sep="")  # NOTE: In this hypothetical example, barcodes 91-96 were not detected
bar.table <- bar.table.full[, good.bars]  # Remove missing bars and summary columns
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))

## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

## Round 2
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

## Round 3
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results3 <- findThresh(call.list=bar.table_sweep.list)
round3.calls <- classifyCells(bar.table, q=findQ(threshold.results3$res, threshold.results3$extrema))
neg.cells <- c(neg.cells, names(round3.calls)[which(round3.calls == "Negative")])

## Repeat until all no negative cells remain (usually 3 rounds)...
final.calls_2 <- c(round3.calls, rep("Negative",length(neg.cells)))
names(final.calls_2) <- c(names(round3.calls),neg.cells)


rm(list = c(ls(pattern = "threshold"),
            ls(pattern = "round"),
            ls(pattern = "bar.table"),
            "neg.cells", "good.bars", "final.calls",
            "i", "n", "q"))


# Rename Barcodes to reflect actual barcode ID
bars <- paste("Bar", 1:11, sep="")
names(bars) <- names(bar.ref)

for (i in 1:length(bars)) {
  final.calls_1[which(final.calls_1 == bars[i])] <- names(bars)[i]
  final.calls_2[which(final.calls_2 == bars[i])] <- names(bars)[i]
}


# Save sample classifications
saveRDS(final.calls_1, "multi_barcodes/final.calls_1.rds")
saveRDS(final.calls_2, "multi_barcodes/final.calls_2.rds")

