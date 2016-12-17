#' @title Print method for echoData Objects.
#' @description Shows main information from echodata Objects.
#'
#' @param x \code{echoData} object provided by \code{readEchograms} function.
#' @param ... Extra argumemts.
#'
#' @export
#' @method print echoData
print.echoData <- function(x, ...){

  cat(paste("\nNumber of echograms: ", x$info$n_echograms, "\n"))

  for(i in seq_along(x$data)){
    cat(paste0("\nFor echogram ", i, ":\n"))

    tempMatrix <- x$data[[i]]

    rangeLon <- c(.getCoordsAxes(min(tempMatrix$dimnames$lon, na.rm = TRUE), "lon"),
                  .getCoordsAxes(max(tempMatrix$dimnames$lon, na.rm = TRUE), "lon"))
    cat(paste("\tRange lon:\tFrom", rangeLon[1], "to", rangeLon[2], "\n"))

    rangeLat <- c(.getCoordsAxes(min(tempMatrix$dimnames$lat, na.rm = TRUE), "lat"),
                  .getCoordsAxes(max(tempMatrix$dimnames$lat, na.rm = TRUE), "lat"))
    cat(paste("\tRange lat:\tFrom", rangeLat[1], "to", rangeLat[2], "\n"))

    rangeTime <- .ac(range(tempMatrix$dimnames$time))
    cat(paste("\tRange time:\tFrom", rangeTime[1], "to", rangeTime[2], "\n"))
  }

  return(invisible())
}

#' @title Summary method for echoData
#' @description Get summary information of echograms included on echodata Objects.
#'
#' @param object \code{echoData} object provided by \code{readEchograms} function.
#' @param ... Extra argumemts.
#'
#' @export
#' @method summary echoData
summary.echoData <- function(object, ...){

  allSummaryData <- list()
  for(i in seq_along(object$data)){

    echoMatriobject <- na.omit(.an(object$data[[i]]$echogram))

    summaryEchogram <- summary(echoMatriobject)

    summaryDimensions <- lapply(object$data[[i]]$dimnames[c("lon", "lat")], summary)
    summaryDimensions <- as.data.frame(do.call(what = "cbind", args = summaryDimensions))
    summaryDimensions <- summaryDimensions[-nrow(summaryDimensions),]

    tempSummary <- data.frame(sA = .an(summaryEchogram),
                              lon = summaryDimensions$lon,
                              lat = summaryDimensions$lat,
                              time = as.POSIXct(as.vector(summary(object$data[[i]]$dimnames$time)),
                                                origin = "1970-01-01 00:00.00 UTC"),
                              stringsAsFactors = FALSE)

    rownames(tempSummary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")

    allSummaryData[[i]] <- tempSummary
  }
  names(allSummaryData) <- paste0("matrix_", seq_along(object$data))

  allSummaryData <- list(n_echograms = object$info$n_echograms,
                         summary_echograms = allSummaryData)

  class(allSummaryData) <- c("summary.echoData")

  return(allSummaryData)
}

#' @title Print method for summary.echoData
#' @description Shows main information from echodata.summary Objects.
#'
#' @param x \code{echoData.summary} object provided by application of summary method to \code{echoData} object.
#' @param ... Extra argumemts.
#'
#' @export
#' @method print summary.echoData
print.summary.echoData <- function(x, ...){

  cat(paste("\nNumber of echograms: ", x$n_echograms, "\n"))

  for(i in seq_along(x$summary_echograms)){
    tempMatrix <- x$summary_echograms[[i]]

    for(j in seq(ncol(tempMatrix))){
      tempMatrix[,j] <- .ac(tempMatrix[,j])
    }

    cat(paste0("\nFor echogram ", i, ":\n"))

    cat(paste(c("\t", colnames(tempMatrix), "\n"), collapse = "\t"))
    for(j in seq(nrow(tempMatrix))){
      cat(paste0("\t", paste(c(c(rownames(tempMatrix)[j], .ac(tempMatrix[j,])), "\n"),
                             collapse = "\t")))
    }
  }

  return(invisible())
}

#' @title Plot method for echoData Objects.
#' @description This function takes an \code{echoData} object and plots an interpolated map
#' showing oxycline depth along shore.
#'
#' @param x \code{echoData} object provided by \code{readEchograms} function.
#' @param ... Arguments passed to \code{\link{echogramPlot}} function.
#'
#' @export
#' @method plot echoData
plot.echoData <- function(x, ...){

  for(i in seq_along(x$data)){

    main <- paste("Echogram", i)

    echogramPlot.default(x = x$data[[i]]$echogram, main = main, ...)

    if(i < length(x$data)){
      readline(prompt = "Hit <enter> to see next echogram:")
    }
  }

  return(invisible())
}

#' @title echogramPlot method for echoData
#'
#' @description This method takes an \code{echoData} object and plots output echograms.
#'
#' @param x Object of class \code{echoData}
#' @param ... Extra arguments passed to \code{\link{echogramPlot}} function.
#'
#' @export
#' @method echogramPlot echoData
echogramPlot.echoData <- function(x, ...){

  for(i in seq_along(x$data)){

    main <- paste("Echogram", i)

    echogramPlot.default(x = x$data[[i]]$echogram, main = main, ...)

    if(i < length(x$data)){
      readline(prompt = "Hit <enter> to see next echogram:")
    }
  }

  return(invisible())
}
