# Filter that removes (converts to NaN) isolated pixels
noiselessFilter <- function(matrixData, radius, times, tolerance){
  radius <- .an(radius)
  times <- .an(times)
  tolerance <- .an(tolerance)

  # Get weighted Matrix (kernel)
  weightedMatrix <- getWeightedMatrix(radius = radius)

  # Apply filters
  finalData <- convolutionQuantile(X = matrixData, kernel = weightedMatrix,
                                   probs = tolerance, times = times)

  return(finalData)
}

# Filter that takes isolated pixels and reforce its closest environment
definerFilter <- function(matrixData, radius, times){
  radius <- .an(radius)
  times <- .an(times)
  tolerance <- 0.1

  # Get weighted Matrix (kernel)
  weightedMatrix <- getWeightedMatrix(radius = radius)

  # Apply filters
  finalData <- convolutionQuantile(X = matrixData, kernel = weightedMatrix,
                                   probs = tolerance, times = times)

  return(finalData)
}

# Function that applies a combination of filters (with different parameters) and
# get a better matrix to calculate limits of oxycline
getFilteredEchogram <- function(fluidMatrix, filterSettings, stepBYstep, ...){

  fluidNames <- dimnames(fluidMatrix$echogram)
  fluidMatrix <- fluidMatrix$echogram

  diffLag <- ifelse(is.null(list(...)$diffLag), 7, list(...)$diffLag)
  diffValues <- apply(fluidMatrix, 2, diff, lag = diffLag)
  diffValues <- abs(unlist(sapply(apply(diffValues, 2, which.max), function(x) if(length(x) < 1) 0 else as.numeric(x)))) #+ diffLag*2

  tempOutput <- sapply(diffValues, function(x, limit) c(rep(1, x), rep(NaN, limit - x)), limit = nrow(fluidMatrix))
  tempOutput <- fluidMatrix*tempOutput

  # Get filtered matrix
  outputList <- list(original = fluidMatrix, tempOutput)
  for(i in seq(nrow(filterSettings))){
    tempFunction <- match.fun(filterSettings[i, "type"])

    tempOutput <- switch(filterSettings[i, "type"],
                         noiselessFilter = tempFunction(matrixData = tempOutput,
                                                         radius = filterSettings[i, "radius"],
                                                         times = filterSettings[i, "times"],
                                                         tolerance = filterSettings[i, "tolerance"]),
                         definerFilter = tempFunction(tempOutput,
                                                       radius = filterSettings[i, "radius"],
                                                       times = filterSettings[i, "times"]),
                         "Incorrect type of filter.")

    dimnames(tempOutput) <- fluidNames

    outputList[[i + 2]] <- tempOutput
  }

  names(outputList) <- c("original", paste0("echogram_", seq(nrow(filterSettings))), "finalEchogram")

  return(outputList)
}

getOxyrange_int <- function(oxyclineData, oxyDims, ...){

  # Define lower and upper limits c(Upper, Lower, Limit)
  if(is.null(list(...)$lineLimits)){
    lineLimits <- c(0.2, 0.99, 0.98)
  }else{
    lineLimits <- list(...)$lineLimits
  }

  allLimits <- list()
  for(i in seq_along(oxyclineData)){
    # Select the final matrix of each echogram and dims
    tempEchogram <- -oxyclineData[[i]]$finalEchogram
    originalMatrix <- -oxyclineData[[i]][[2]]

    # Replace NA with zeros
    tempEchogram[is.na(tempEchogram)] <- 0
    originalMatrix[is.na(originalMatrix)] <- 0

    # Extract original values from filtered matrix
    cumsumMatrix <- apply(tempEchogram, 2, function(x) cumsum(x)/sum(x))
    cumsumMatrix <- apply(cumsumMatrix, 2, findInterval, vec = c(-Inf, lineLimits[1:2], Inf))
    tempEchogram[cumsumMatrix != 2] <- 0

    # Extract oxycline limits
    cumsumMatrix <- apply(tempEchogram, 2, function(x) cumsum(x)/sum(x))
    cumsumMatrix <- apply(cumsumMatrix, 2, findInterval, vec = c(-Inf, lineLimits[3], Inf))

    limitsData <- rep(NA, ncol(cumsumMatrix))

    index <- colSums(cumsumMatrix, na.rm = TRUE) > 0
    limitsData[index] <- apply(cumsumMatrix[,index], 2, function(x) na.omit(which(!duplicated(x))))[2,]
    limitsData[index] <- sapply(limitsData[index], function(x, depths) depths[x],
                                depths = .an(rownames(tempEchogram)))

    # Make smoothing
    limitsData <- smoothVector(limitsData, spar = ifelse(is.null(list(...)$spar), 0.1, list(...)$spar))
    # limitsData <- movingAverage(x = limitsData, n = 20)

    # Reorganize info
    limitsData <- data.frame(lon = oxyDims[[i]]$lon,
                             lat = oxyDims[[i]]$lat,
                             oxycline_limit = limitsData,
                             stringsAsFactors = FALSE)
    rownames(limitsData) <- colnames(tempEchogram)

    # Compile values on a list
    allLimits[[i]] <- limitsData
  }

  names(allLimits) <- names(oxyclineData)

  return(allLimits)
}

# Get dimensions from echogram data
getOxyDims <- function(oxyclineData){
  allDims <- list()
  for(i in seq_along(oxyclineData)){
    tempDim <- list(lon = oxyclineData[[i]]$dimnames$lon,
                    lat = oxyclineData[[i]]$dimnames$lat,
                    time = oxyclineData[[i]]$dimnames$time)

    allDims[[i]] <- tempDim
  }

  names(allDims) <- paste0("matrix_", seq_along(oxyclineData))

  return(allDims)
}

getEchoData <- function(fileMode, directoryMode,
                         validFish38, validBlue38, upLimitFluid120,
                         pinInterval, date.format, ...){

  # dayHours <- c("06:30", "17:30")
  # nightHours <- c("19:30", "04:50")

  if(is.null(fileMode) & is.null(directoryMode)){
    stop("At least whether fileMode or directoryMode must not be NULL.")
  }

  if(!is.null(fileMode)){

    # Check if files exist or not
    if(all(!sapply(fileMode, file.exists))){
      stop("Some given files do not exist.")
    }

    # Read Matlab files (.m)
    fish38_matrix <- readMat(fileMode$fish38_file)
    fluid120_matrix <- readMat(fileMode$fluid120_file)
    blue38_matrix <- readMat(fileMode$blue38_file)

    # Extract lon, lat, time and depth
    allLon <- fluid120_matrix$Longitude
    allLat <- fluid120_matrix$Latitude
    allTime <- paste(.ac(fluid120_matrix$Ping.date),
                     .ac(fluid120_matrix$Ping.time))
    depth <- fluid120_matrix$depth

    # Extract echogram matrix
    fish38_matrix <- fish38_matrix$Data.values
    fluid120_matrix <- fluid120_matrix$Data.values
    blue38_matrix <- blue38_matrix$Data.values

    # Clear data using limit parameters
    index <- (fish38_matrix < -998 | fish38_matrix < validFish38[1] |
                fish38_matrix > validFish38[2])
    fish38_matrix[index] <- NaN

    index <- (blue38_matrix < -998 | blue38_matrix < validBlue38[1] |
                blue38_matrix > validBlue38[2])
    blue38_matrix[index] <- NaN

    index <- fluid120_matrix < -998 | fluid120_matrix > upLimitFluid120
    fluid120_matrix[index] <- NaN

    # # Clear main data (Fluid-like) using Fish and Blue noise data
    # allData <- fluid120_matrix*(is.na(blue38_matrix) & is.na(fish38_matrix))
    # allData[allData == 0] <- NaN

    # Make operations between matrices
    if(is.null(list(...)$operationType)){
      operationType <- 1
    }else{
      operationType <- list(...)$operationType
    }

    # operationType <- 1
    operationsInput <- list(fluid120_matrix = fluid120_matrix,
                            blue38_matrix = blue38_matrix,
                            fish38_matrix = fish38_matrix)

    # Apply operations between matrices
    allData <- operationsSet(operationsInput = operationsInput, type = operationType)

    # Transpose matrix
    allData <- t(allData)

  }else{

    directory <- directoryMode$directory

    # Define text pattern of data bases
    pattern_Fish38  <- directoryMode$fish38_pattern
    pattern_Blue38  <- directoryMode$blue38_pattern
    pattern_Fluid120  <- directoryMode$fluid120_pattern

    # Generate file list with text patterns
    listFiles_Fish <- list.files(path = directory, pattern = pattern_Fish38,
                                 full.names = TRUE, recursive = TRUE)
    listFiles_Fluid <- list.files(path = directory, pattern = pattern_Fluid120,
                                  full.names = TRUE, recursive = TRUE)
    listFiles_Blue <- list.files(path = directory, pattern = pattern_Blue38,
                                 full.names = TRUE, recursive = TRUE)

    # Read files and concatenate in one matrix
    allData <- allTime <- allLon <- allLat <- NULL
    for(i in seq_along(listFiles_Fish)){
      tempList_Fish <- readMat(listFiles_Fish[i])
      tempList_Fluid <- readMat(listFiles_Fluid[i])
      tempList_Blue <- readMat(listFiles_Blue[i])

      tempData_Fish <- tempList_Fish$Data.values
      tempData_Fluid <- tempList_Fluid$Data.values
      tempData_Blue <- tempList_Blue$Data.values

      if(i == 1)
        depth <- .an(tempList_Fluid$depth)

      tempTime <- paste(.ac(tempList_Fluid$Ping.date),
                        .ac(tempList_Fluid$Ping.time))
      tempLon <- tempList_Fluid$Longitude
      tempLat <- tempList_Fluid$Latitude
      rm(list = c("tempList_Fish", "tempList_Fluid", "tempList_Blue"))

      # Clear data using limit parameters
      tempData_Fish[tempData_Fish < -998 | tempData_Fish < validFish38[1] |
                      tempData_Fish > validFish38[2]] <- NaN
      tempData_Blue[tempData_Blue < -998 | tempData_Blue < validBlue38[1] |
                      tempData_Blue > validBlue38[2]] <- NaN
      tempData_Fluid[tempData_Fluid < -998 | tempData_Fluid > upLimitFluid120] <- NaN

      # Clear main data (Fluid-like) using Fish and Blue noise data
      # tempData <- tempData_Fluid*(is.na(tempData_Blue) & is.na(tempData_Fish))
      # tempData[tempData == 0] <- NaN

      # Make operations between matrices
      if(is.null(list(...)$operationType)){
        operationType <- 1
      }else{
        operationType <- list(...)$operationType
      }

      # operationType <- 1
      operationsInput <- list(fluid120_matrix = tempData_Fluid,
                              blue38_matrix = tempData_Blue,
                              fish38_matrix = tempData_Fish)

      # Apply operations between matrices
      tempData <- operationsSet(operationsInput = operationsInput, type = operationType)

      allLon <- c(allLon, tempLon)
      allLat <- c(allLat, tempLat)
      allTime <- c(allTime, tempTime)
      allData <- cbind(allData, t(tempData))
    }
  }

  # Convert time
  allTime <- strptime(allTime, format = date.format)

  if(sum(is.na(allTime)) > 0)
    stop("Incorrect value for 'date.format'.")

  # Get points where the difference between two pin is larger than pinInterval (sec)
  breakPoints <- which(.an(diff(allTime)) > pinInterval)

  index <- is.null(dim(breakPoints)) && length(breakPoints) > 1
  breakPoints <- c(0, if(index) NULL else breakPoints, dim(allData)[2])

  # Split big matrix by breakpoints to get matrix of echograms
  data <- list()
  for(i in seq(2, length(breakPoints))){
    tempEchogram <- list()

    # Get index by breakpoints
    index <- seq(breakPoints[i - 1] + 1, breakPoints[i])

    # Split data by index
    tempMatrix <- allData[,index]
    tempTimes <- allTime[index]
    tempLon <- allLon[index]
    tempLat <- allLat[index]

    colnames(tempMatrix) <- .ac(tempTimes)
    rownames(tempMatrix) <- round(depth, 2)

    tempEchogram[[1]] <- tempMatrix
    tempEchogram[[2]] <- list(depth = depth,
                              time = tempTimes,
                              lon = tempLon,
                              lat = tempLat)

    names(tempEchogram) <- c("echogram", "dimnames")

    data[[i - 1]] <- tempEchogram
  }

  names(data) <- paste0("matrix_", seq_along(breakPoints[-1]))

  # Build final output
  output <- list(info = list(parameters = list(validFish38 = validFish38,
                                               validBlue38 = validBlue38,
                                               upLimitFluid120 = upLimitFluid120),
                             n_echograms = length(breakPoints) - 1),
                 data = data)

  return(output)
}

.echogramPlot <- function(echogram, colEchogram, cex.axis = 1, ...){

  # Define raster from echogram
  xAxis <- as.POSIXct(dimnames(echogram)[[2]])
  yAxis <- abs(.an(dimnames(echogram)[[1]]))

  ext_xAxis <- seq.POSIXt(range(xAxis)[1], range(xAxis)[2], by = "sec")

  newEchogram <- matrix(NA, nrow = nrow(echogram), ncol = length(ext_xAxis))
  newEchogram[,match(as.integer(xAxis), as.integer(ext_xAxis))] <- echogram

  newEchogram <- t(newEchogram)
  newEchogram <- newEchogram[,ncol(newEchogram):1]

  # Get plot of raster
  nIntervals <- 5

  xVector <- seq.POSIXt(from = min(xAxis), to = max(xAxis), length.out = nIntervals)

  xlim <- range(xVector)
  ylim <- range(pretty(yAxis, n = nIntervals))

  par(mar = c(3, 4, 2, 3), xaxs = "i", yaxs = "i")

  image(x = ext_xAxis, y = yAxis, z = newEchogram,
        xlim = .an(xlim), ylim = ylim, axes = FALSE, ylab = "Depth (m)",
        useRaster = FALSE, col = colEchogram, ...)

  axis(2, at = pretty(yAxis), labels = rev(abs(pretty(yAxis))), las = 2, cex.axis = cex.axis)
  axis(1, at = .an(xVector), labels = as.Date(xVector), cex.axis = cex.axis)
  axis(1, at = .an(xVector), labels = strftime(xVector, format="%H:%M:%S"), line = 1,
       tick = FALSE, cex.axis = cex.axis)

  box()

  return(invisible())
}
