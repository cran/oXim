#' @title Print method for oxyclineData Objects.
#' @description Shows main information from oxyclineData Objects.
#'
#' @param x Object of class \code{oxyclineData}.
#' @param ... Extra argumemts.
#'
#' @export
#' @method print oxyclineData
print.oxyclineData <- function(x, ...){

  for(i in seq_along(x$oxycline_range)){

    tempData <- x$oxycline_range[[i]]

    cat(paste0("\nFor echogram ", i, ":\n"))

    rangeLon <- c(.getCoordsAxes(min(tempData$lon, na.rm = TRUE), "lon"),
                  .getCoordsAxes(max(tempData$lon, na.rm = TRUE), "lon"))
    cat(paste("\tRange lon:\tFrom", rangeLon[1], "to", rangeLon[2], "\n"))

    rangeLat <- c(.getCoordsAxes(min(tempData$lat, na.rm = TRUE), "lat"),
                  .getCoordsAxes(max(tempData$lat, na.rm = TRUE), "lat"))
    cat(paste("\tRange lat:\tFrom", rangeLat[1], "to", rangeLat[2], "\n"))

    rangeTime <- .ac(range(as.POSIXct(rownames(tempData))))
    cat(paste("\tRange time:\tFrom", rangeTime[1], "to", rangeTime[2], "\n"))

    rangeOxylimits <- round(abs(range(tempData$oxycline_limit, na.rm = TRUE)), 1)
    cat(paste("\tRange oxycline depth:\tFrom", rangeOxylimits[2], "m to", rangeOxylimits[1], "m\n"))
  }

  return(invisible())
}

#' @title Summary method for oxyclineData Objects.
#' @description Get summary information of oxycline depth limits included on \code{oxyclineData} objects.
#'
#' @param object Object of class \code{oxyclineData}.
#' @param ... Extra argumemts.
#'
#' @export
#' @method summary oxyclineData
summary.oxyclineData <- function(object, ...){

  allSummaryData <- list()
  for(i in seq_along(object$oxycline_range)){

    tempData <- object$oxycline_range[[i]]

    # Get summary for oxycline depth info
    oxylimitsObject <- na.omit(tempData$oxycline_limit)
    summaryOxylimits <- .an(summary(oxylimitsObject))

    # Get summary for coords info
    summaryCoords <- data.frame(lon = .an(summary(na.omit(tempData$lon))),
                                lat = .an(summary(na.omit(tempData$lat))),
                                stringsAsFactors = FALSE)

    # Get summary for time info
    summaryTime <- .ac(summary(as.POSIXct(rownames(tempData))))

    # Join all summary info in one single data frame
    tempSummary <- data.frame(lon = summaryCoords$lon, lat = summaryCoords$lat,
                              limits = summaryOxylimits, time = summaryTime,
                              stringsAsFactors = FALSE)

    rownames(tempSummary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")

    allSummaryData[[i]] <- tempSummary
  }
  names(allSummaryData) <- paste0("matrix_", seq_along(object$data))

  allSummaryData <- list(n_database = length(object$oxycline_range),
                         summary_database = allSummaryData)

  class(allSummaryData) <- c("summary.oxyclineData")

  return(allSummaryData)
}

#' @title Print method for summary.oxyclineData
#' @description Shows main information from \code{oxyclineData.summary} objects.
#'
#' @param x Object of class \code{summary.oxyclineData}.
#' @param ... Extra argumemts.
#'
#' @export
#' @method print summary.oxyclineData
print.summary.oxyclineData <- function(x, ...){

  cat(paste("\nNumber of database: ", x$n_database, "\n"))

  for(i in seq_along(x$summary_database)){
    tempMatrix <- x$summary_database[[i]]

    for(j in seq(ncol(tempMatrix))){
      tempMatrix[,j] <- .ac(tempMatrix[,j])
    }

    cat(paste0("\nFor database ", i, ":\n"))

    cat(paste(c("\t", colnames(tempMatrix), "\n"), collapse = "\t"))
    for(j in seq(nrow(tempMatrix))){
      cat(paste0("\t", paste(c(c(rownames(tempMatrix)[j], .ac(tempMatrix[j,])), "\n"),
                             collapse = "\t")))
    }
  }

  return(invisible())
}

#' @title Plot method for oxyclineData
#'
#' @description This method takes an \code{oxyclineData} object, make an interpolation of oxycline
#' values and show them on a map. Interpolation methodology is based on \code{akima} package.
#'
#' @param x Object of class \code{oxyclineData}
#' @param interpParams \code{list} object including parameters passed to \code{\link{idw}} function.
#' @param xlengthAxes Desired length of the axis 'x' labels.
#' @param ylengthAxes Desired length of the axis 'y' labels.
#' @param ... Extra arguments passed to \code{\link{image}} function.
#' @param showLimits \code{logical}. Do you want to show the oxycline limit line?
#'
#' @export
#' @method plot oxyclineData
plot.oxyclineData <- function(x, interpParams = list(myGrid = NULL), xlengthAxes = 6, ylengthAxes = 6,
                              showLimits = FALSE, ...){

  # Combine all matrices in one data.frame
  allData <- NULL
  for(i in seq_along(x$oxycline_range)){
    allData <- rbind(allData, x$oxycline_range[[i]])
  }

  # Change colnames
  allData <- as.data.frame(allData, stringsAsFactors = FALSE)
  colnames(allData) <- colnames(x$oxycline_range[[1]])
  x <- allData

  # Remove rows with no-data in lat, lon or upper_limit
  index <- complete.cases(x[,c("lon", "lat", "oxycline_limit")])
  x <- x[index,]

  # Make interpolation using akima package
  x <- do.call(what = ".interpIDW",
               args = c(list(myData = x, XYZnames = c("lon", "lat", "oxycline_limit")), interpParams))

  myGrid <- x$myGrid
  x <- x$myIDW

  x <- list(x = unique(x@coords[,1]), y = unique(x@coords[,2]),
            z = matrix(data = x@data$var1.pred, nrow = x@grid@cells.dim[1], ncol = x@grid@cells.dim[2]))

  # Plot map
  originalPar <- par()

  par(mar = c(4, 4, 1, 1))
  image(x = x, axes = FALSE, ...)

  xAxes <- range(x$x)
  xAxes <- seq(from = xAxes[1], to = xAxes[2], length.out = xlengthAxes)

  yAxes <- range(x$y)
  yAxes <- seq(from = yAxes[1], to = yAxes[2], length.out = ylengthAxes)

  if(!is.null(list(...)$cex.axis)){
    axis(side = 1, at = xAxes, labels = sapply(xAxes, .getCoordsAxes, what = "lon"),
         cex.axis = list(...)$cex.axis)
    axis(side = 2, at = yAxes, labels = sapply(yAxes, .getCoordsAxes, what = "lat"),
         las = 2, cex.axis = list(...)$cex.axis)
  }else{
    axis(side = 1, at = xAxes, labels = sapply(xAxes, .getCoordsAxes, what = "lon"))
    axis(side = 2, at = yAxes, labels = sapply(yAxes, .getCoordsAxes, what = "lat"), las = 2)
  }

  box()

  return(invisible())
}

#' @title echogramPlot method for oxyclineData
#'
#' @description This method takes an \code{oxyclineData} object and plots output echograms. Optionaly,
#' users can add oxycline line.
#'
#' @param x Object of class \code{oxyclineData} with internal echogram matrix to be plotted.
#' @param colEchogram Pallete of colours to plot the echograms. If \code{NULL} (default) the system
#' will use the same combination used on object \code{colPallete}.
#' @param oxyLine \code{logical}. Do you want to add oxycline line to the plot?
#' @param oxyLineParams If \code{oxyLine = TRUE}, parameters passed to \code{\link{lines}} function.
#' @param ... Extra arguments passed to \code{\link{echogramPlot}} function.
#'
#' @method echogramPlot oxyclineData
#' @export
echogramPlot.oxyclineData <- function(x, colEchogram = "colPalette", oxyLine = TRUE, oxyLineParams = list(), ...){

  for(i in seq_along(x$outputs)){
    # Plot echogram
    echogramPlot.default(x = x$outputs[[i]]$original, colEchogram = colEchogram, ...)

    # Add oxycline line
    if(isTRUE(oxyLine)){

      timeVector <- as.POSIXct(rownames(x$oxycline_range[[i]]))
      oxyLimit <- x$oxycline_range[[i]]$oxycline_limit
      oxyLimit <- abs(min(.an(rownames(x$outputs[[i]]$finalEchogram)))) + oxyLimit

      do.call(what = "lines",
              args = c(with(x$oxycline_range[[i]], list(x = timeVector, y = oxyLimit)), oxyLineParams))
    }
  }

  return(invisible())
}


#' @rdname echogramPlot
#' @method echogramPlot matrix
#' @export
echogramPlot.matrix <- function(x, ...){

  echogramPlot.default(x = x, ...)

  return(invisible())
}
