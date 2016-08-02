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

    rangeOxylimits <- round(abs(range(tempData$upper_limit, na.rm = TRUE)), 1)
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
    oxylimitsObject <- na.omit(tempData$upper_limit)
    summaryOxylimits <- .an(summary(oxylimitsObject))

    # Get summary for coords info
    summaryCoords <- data.frame(lon = .an(summary(tempData$lon)),
                                lat = .an(summary(tempData$lat)),
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
#'
#' @export
#' @method plot oxyclineData
plot.oxyclineData <- function(x, interpParams = list(myGrid = NULL), xlengthAxes = 5, ylengthAxes = 5, ...){

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
  index <- complete.cases(x[,c("lon", "lat", "upper_limit")])
  x <- x[index,]

  # Make interpolation using akima package
  x <- do.call(what = ".interpIDW",
               args = c(list(myData = x, XYZnames = c("lon", "lat", "upper_limit")), interpParams))

  myGrid <- x$myGrid
  x <- x$myIDW

  x <- matrix(data = x@data$var1.pred, nrow = x@grid@cells.dim[1], ncol = x@grid@cells.dim[2])

  # Plot map
  originalPar <- par()

  par(mar = c(4, 4, 1, 1))
  image(x = x, axes = FALSE, ...)

  xAxes <- range(myGrid$x)
  xAxes <- seq(from = xAxes[1], to = xAxes[2], length.out = xlengthAxes)

  yAxes <- range(myGrid$y)
  yAxes <- seq(from = yAxes[1], to = yAxes[2], length.out = ylengthAxes)

  axis(side = 1, at = seq(from = 0, to = 1, length.out = xlengthAxes),
       labels = sapply(xAxes, .getCoordsAxes, what = "lon"))
  axis(side = 2, at = seq(from = 0, to = 1, length.out = ylengthAxes),
       labels = sapply(yAxes, .getCoordsAxes, what = "lat"), las = 2)
  box()

  return(invisible())
}
