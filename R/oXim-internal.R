#' @rdname echogramPlot
#' @export
echogramPlot.default <- function(x, colEchogram = "colPalette", ...){
  colEchogram <- get(colEchogram)
  .echogramPlot(echogram = x, colEchogram = colEchogram, ...)

  return(invisible())
}

operationsSet <- function(operationsInput, type){

  operationFunction <- match.fun(paste0("operationFunction_", type))

  output <- operationFunction(operationsInput)

  return(output)
}

operationFunction_1 <- function(operationsInput){

  # Check input data
  neededMatrices <- c("fluid120_matrix", "blue38_matrix", "fish38_matrix")

  if(is.list(operationsInput) && all(is.element(neededMatrices, names(operationsInput)))){

    # Check classes
    allClasses <- unique(sapply(operationsInput, class))

    if(length(allClasses) > 1 || allClasses != "matrix"){
      stop("'operationsInput' must be a list of echogram matrices.")
    }

    # Check if some matrix have different dimensions
    index <- apply(sapply(operationsInput, dim), 1, function(x) length(unique(x)))
    if(any(index != 1)){
      stop("'operationsInput' must be a list of echogram matrices with the same dimensions.")
    }

  }else{
    stop("'operationsInput' must be a list of echogram matrices.")
  }

  # Make operations
  for(i in seq_along(operationsInput)){
    operationsInput[[i]][is.na(operationsInput[[i]])] <- 0
  }

  output <- with(operationsInput, ((fluid120_matrix*1e3) + (blue38_matrix*1e3) + fish38_matrix)/3)

  return(output)
}
