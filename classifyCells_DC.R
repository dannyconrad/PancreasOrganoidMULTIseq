classifyCells_DC <- function (barTable, q) 
{
  require(KernSmooth)
  barTable.n <- as.data.frame(log2(barTable))
  for (i in 1:ncol(barTable)) {
    ind <- which(is.finite(barTable.n[, i]) == FALSE)
    barTable.n[ind, i] <- 0
    barTable.n[, i] <- barTable.n[, i] - mean(barTable.n[, 
                                                         i])
  }
  n_BC <- ncol(barTable.n)
  n_cells <- nrow(barTable.n)
  bc_calls <- rep("Negative", n_cells)
  names(bc_calls) <- rownames(barTable.n)
  for (i in 1:n_BC) {
    thresh <- NULL
    model <- tryCatch({
      approxfun(bkde(barTable.n[, i], kernel = "normal"))
    }, error = function(e) {
      print(paste0("No threshold found for ", colnames(barTable.n)[i], 
                   "..."))
    })
    if (class(model) == "character") {
      next
    }
    x <- seq(from = quantile(barTable.n[, i], 0.001), to = quantile(barTable.n[, 
                                                                               i], 0.999), length.out = 100)
    extrema <- localMaxima(model(x))
    if (length(extrema) <= 1) {
      print(paste0("No threshold found for ", colnames(barTable.n)[i], 
                   "..."))
      hist(barTable.n[,i], breaks = 100)
      thresh <- as.numeric(readline(prompt="Please manually choose a threshold: "))
    }

    if (is.null(thresh)) {
      low.extreme <- extrema[which.max(model(x)[extrema])]
      high.extreme <- max(extrema)
      if (low.extreme == high.extreme) {
        print(paste0("No threshold found for ", colnames(barTable.n)[i], 
                     "..."))
      }
      thresh <- quantile(c(x[high.extreme], x[low.extreme]), 
                         q)     
    }
    
    cell_i <- which(barTable.n[, i] >= thresh)
    n <- length(cell_i)
    if (n == 0) {
      next
    }
    bc_calls[cell_i] <- sapply(bc_calls[cell_i], FUN = function(x) {
      if (x == "Negative") {
        return(colnames(barTable.n)[i])
      }
      else {
        return("Doublet")
      }
    })
  }
  return(bc_calls)
}