# ordfiltInR <- function(data, x, weightedMatrix){
#   newData <- data
#
#   # Corner UpLeft
#   i <- 1
#   j <- 1
#   miniData <- data[(i):(i + 2), (j):(j + 2)] * weightedMatrix
#   newData[i, j] <- miniData[order(miniData)[x]]
#
#   # Corner DownLeft
#   i <- nrow(data)
#   j <- 1
#   miniData <- data[(i):(i - 2), (j):(j + 2)] * weightedMatrix
#   newData[i, j] <- miniData[order(miniData)[x]]
#
#   # Corner UpRight
#   i <- 1
#   j <- ncol(data)
#   miniData <- data[(i):(i + 2), (j):(j - 2)] * weightedMatrix
#   newData[i, j] <- miniData[order(miniData)[x]]
#
#   # Corner DownRight
#   i <- nrow(data)
#   j <- ncol(data)
#   miniData <- data[(i):(i - 2), (j):(j - 2)] * weightedMatrix
#   newData[i, j] <- miniData[order(miniData)[x]]
#
#
#   # Downside
#   i <- nrow(data):(nrow(data) - 2)
#   for(j in seq(from = 2, to = ncol(data) - 1)){
#     miniData <- data[i, (j - 1):(j + 1)] * weightedMatrix
#     newData[i, (j - 1):(j + 1)] <- miniData[order(miniData)[x]]
#   }
#
#   # Leftside
#   j <- 1:3
#   for(i in seq(from = 2, to = nrow(data) - 1)){
#     miniData <- data[(i - 1):(i + 1), j] * weightedMatrix
#     newData[(i - 1):(i + 1), j] <- miniData[order(miniData)[x]]
#   }
#
#   # Upside
#   i <- 1:3
#   for(j in seq(from = 2, to = ncol(data) - 1)){
#     miniData <- data[i, (j - 1):(j + 1)] * weightedMatrix
#     newData[i, (j - 1):(j + 1)] <- miniData[order(miniData)[x]]
#   }
#
#   # Rightside
#   j <- ncol(data):(ncol(data) - 2)
#   for(i in seq(from = 2, to = nrow(data) - 1)){
#     miniData <- data[(i - 1):(i + 1), j] * weightedMatrix
#     newData[(i - 1):(i + 1), j] <- miniData[order(miniData)[x]]
#   }
#
#   # No borders
#   for(j in seq(from = 2, to = ncol(data) - 1)){
#     for(i in seq(from = 2, to = nrow(data) - 1)){
#       miniData <- data[(i - 1):(i + 1), (j - 1):(j + 1)]
#
#       if(sum(!is.na(miniData), na.rm = TRUE) == 0)
#         next
#
#       newData[i, j] <- miniData[order(miniData * weightedMatrix)[x]]
#     }
#   }
#
#   return(newData)
# }

.ordfilt2_C <- function(data, x, weightedMatrix){

  # Error messages
  if(mode(data) != "numeric" | mode(weightedMatrix) != "numeric")
    stop("Incorrect mode of data or weightedMatrix (both must be 'numeric').")

  # No borders
  miniData <- ordfiltInC(data = data, x = as.integer(x), weightedMatrix = .an(weightedMatrix))

  return(miniData)
}

# Filter that removes (converts to NaN) isolated pixels
.noiselessFilter <- function(data, radius, times, tolerance){
  radius <- .an(radius)
  times <- .an(times)
  tolerance <- .an(tolerance)

  # Get range of values
  rangeValues <- range(data, na.rm = TRUE)

  # Convert NA to -999
  data[is.na(data)] <- -999

  # Get weighted Matrix
  weightedMatrix <- diag(radius) + diag(radius)[,radius:1]
  constant1 <- ceiling(radius/2)
  weightedMatrix[constant1,] <- 1

  constant2 <- ceiling(sum(weightedMatrix)*tolerance)

  if(constant2 - 1 < 0)
    stop("Incorrect value for tolerance.")

  finalData <- data
  for(i in 1:times)
    finalData <- .ordfilt2_C(data = finalData,
                             x = constant2,
                             weightedMatrix = weightedMatrix)

  finalData[finalData < rangeValues[1] | finalData == 0] <- NA

  return(finalData)
}

# Filter that takes isolated pixels and reforce its closest environment
.definerFilter <- function(data, radius, times){
  radius <- .an(radius)
  times <- .an(times)

  # Get range of values
  rangeValues <- range(data, na.rm = TRUE)

  # Convert NA to 999
  data[is.na(data)] <- 999

  # Get weighted Matrix
  weightedMatrix <- diag(radius) + diag(radius)[,radius:1]
  constant1 <- ceiling(radius/2)
  weightedMatrix[constant1,] <- 1

  constant2 <- 1

  finalData <- data
  for(i in 1:times)
    finalData <- .ordfilt2_C(data = finalData,
                             x = constant2,
                             weightedMatrix = weightedMatrix)

  finalData[finalData > rangeValues[2]] <- NA

  return(finalData)
}

# Function that applies a combination of filters (with different parameters) and
# get a better matrix to calculate limits of oxycline
.getFilteredEchogram <- function(fluidMatrix, filterSettings, stepBYstep){

  fluidNames <- dimnames(fluidMatrix$echogram)
  fluidMatrix <- fluidMatrix$echogram

  # Get filtered matrix
  tempOutput <- fluidMatrix
  outputList <- list(original = fluidMatrix)
  for(i in seq(nrow(filterSettings))){
    tempFunction <- match.fun(filterSettings[i, "type"])

    tempOutput <- switch(filterSettings[i, "type"],
                         .noiselessFilter = tempFunction(tempOutput,
                                                         radius = filterSettings[i, "radius"],
                                                         times = filterSettings[i, "times"],
                                                         tolerance = filterSettings[i, "tolerance"]),
                         .definerFilter = tempFunction(tempOutput,
                                                       radius = filterSettings[i, "radius"],
                                                       times = filterSettings[i, "times"]),
                         "Incorrect type of filter.")

    dimnames(tempOutput) <- fluidNames

    if(stepBYstep)
      outputList[[i + 1]] <- tempOutput else
        if(i == nrow(filterSettings))
          outputList[[2]] <- tempOutput
  }

  names(outputList) <- if(length(outputList) > 2)
    c("original", paste0("echogram_", seq(nrow(filterSettings) - 1)), "finalEchogram") else
      c("original", "finalEchogram")

  return(outputList)
}

.getOxyrange <- function(oxyclineData, oxyDims){

  allLimits <- list()
  for(i in seq_along(oxyclineData)){
    # Select the final matrix of each echogram and dims
    tempEchogram <- oxyclineData[[i]]$finalEchogram
    tempDims <- oxyDims[[i]]

    # Define lower and upper limits
    lineLimits <- c(0.10, 0.98)

    # Get matrix where values of tempEchogram are lower than zero
    dataEchogram <- drop(outer(tempEchogram, 0, "<"))
    dataEchogram[is.na(dataEchogram)] <- 0

    # What columns has, at least, one value for getting oxycline range
    index <- which(colSums(dataEchogram) > 0)

    # Get sums by column
    sumByCol <- apply(dataEchogram, 2, sum)

    lineLimits <- ceiling(sapply(lineLimits, "*", sumByCol))

    dataEchogram <- apply(dataEchogram, 2, cumsum)

    # Set empty matrix for recording range values
    limitsData <- matrix(NA, nrow = ncol(tempEchogram), ncol = 4)
    dimnames(limitsData) <- list(colnames(tempEchogram),
                                 c("lower_limit", "upper_limit", "lon", "lat"))
    for(j in index){
      # Get vector with cumsum values by depth
      limitIndex <- !duplicated(dataEchogram[,j])
      limitVector <- dataEchogram[limitIndex, j]

      # Select and save limit values
      limitIndex <- match(as.numeric(lineLimits[j,]), limitVector)
      limitsData[j, c(1, 2)] <- .an(names(limitVector)[limitIndex])
      # c(, tempDims$lon[j], tempDims$lat[j])
    }

    limitsData[,c(3, 4)] <- cbind(tempDims$lon, tempDims$lat)

    # Compile values on a list
    output <- as.data.frame(limitsData, stringsAsFactors = FALSE)
    allLimits[[i]] <- output[complete.cases(output[,c("lon", "lat")]),]
  }

  names(allLimits) <- names(oxyclineData)

  return(allLimits)
}

# Get dimensions from echogram data
.getOxyDims <- function(oxyclineData){
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

# Takes outputs from Echopen and generates a matrix to calculate Oxycline
# .getEchoData <- function(directory, validFish38, validBlue38, upLimitFluid120, pinInterval,
#                          date.format){
#
#   # Define ttext pattern of databases
#   pattern_Fish38  <- "_2Freq_Fish38.mat"
#   # pattern_Fluid38 <- "_2Freq_Fluid38.mat"
#   pattern_Blue38  <- "_2Freq_Blue38.mat"
#
#   # pattern_Fish120   <- "_2Freq_Fish120.mat"
#   pattern_Fluid120  <- "_2Freq_Fluid120.mat"
#   # pattern_Blue120   <- "_2Freq_Blue120.mat"
#
#   # Generate file list with text patterns
#   listFiles_Fish <- list.files(path = directory, pattern = pattern_Fish38,
#                                full.names = TRUE, recursive = TRUE)
#   listFiles_Fluid <- list.files(path = directory, pattern = pattern_Fluid120,
#                                 full.names = TRUE, recursive = TRUE)
#   listFiles_Blue <- list.files(path = directory, pattern = pattern_Blue38,
#                                full.names = TRUE, recursive = TRUE)
#
#   # Read files and concatenate in one matrix
#   allData <- allTime <- allLon <- allLat <- NULL
#   for(i in seq_along(listFiles_Fish)){
#     tempList_Fish <- readMat(listFiles_Fish[i])
#     tempList_Fluid <- readMat(listFiles_Fluid[i])
#     tempList_Blue <- readMat(listFiles_Blue[i])
#
#     tempData_Fish <- tempList_Fish$Data.values
#     tempData_Fluid <- tempList_Fluid$Data.values
#     tempData_Blue <- tempList_Blue$Data.values
#
#     if(i == 1)
#       depth <- as.numeric(tempList_Fluid$depth)
#
#     tempTime <- paste(as.character(tempList_Fluid$Ping.date),
#                       as.character(tempList_Fluid$Ping.time))
#     tempLon <- tempList_Fluid$Longitude
#     tempLat <- tempList_Fluid$Latitude
#     rm(list = c("tempList_Fish", "tempList_Fluid", "tempList_Blue"))
#
#     # Clear data using limit parameters
#     tempData_Fish[tempData_Fish < -998 | tempData_Fish < validFish38[1] |
#                     tempData_Fish > validFish38[2]] <- NaN
#     tempData_Blue[tempData_Blue < -998 | tempData_Blue < validBlue38[1] |
#                     tempData_Blue > validBlue38[2]] <- NaN
#     tempData_Fluid[tempData_Fluid < -998 | tempData_Fluid > upLimitFluid120] <- NaN
#
#     # Clear main data (Fluid-like) using Fish and Blue noise data
#     tempData <- tempData_Fluid*(is.na(tempData_Blue) & is.na(tempData_Fish))
#     tempData[tempData == 0] <- NaN
#
#     allLon <- c(allLon, tempLon)
#     allLat <- c(allLat, tempLat)
#     allTime <- c(allTime, tempTime)
#     allData <- cbind(allData, t(tempData))
#   }
#
#   # Convert time
#   allTime <- strptime(allTime, format = date.format)
#
#   if(sum(is.na(allTime)) > 0)
#     stop("Incorrect value for 'date.format'.")
#
#   # Get points where the difference between two pin is larger than pinInterval (sec)
#   breakPoints <- which(as.numeric(diff(allTime)) > pinInterval)
#   breakPoints <- if(is.null(dim(breakPoints)) && length(breakPoints) > 1)
#     c(0, dim(allData)[2]) else
#       c(0, breakPoints, dim(allData)[2])
#
#   # Split big matrix by breakpoints to get matrix of echograms
#   data <- list()
#   for(i in seq(2, length(breakPoints))){
#     tempEchogram <- list()
#
#     index <- seq(breakPoints[i - 1] + 1, breakPoints[i])
#
#     tempMatrix <- allData[,index]
#     tempTimes <- allTime[index]
#     tempLon <- allLon[index]
#     tempLat <- allLat[index]
#
#     colnames(tempMatrix) <- as.character(tempTimes)
#     rownames(tempMatrix) <- round(depth, 2)
#
#     tempEchogram[[1]] <- tempMatrix
#     tempEchogram[[2]] <- list(depth = depth,
#                               time = tempTimes,
#                               lon = tempLon,
#                               lat = tempLat)
#
#     names(tempEchogram) <- c("echogram", "dimnames")
#
#     data[[i - 1]] <- tempEchogram
#   }
#
#   names(data) <- paste0("matrix_", seq_along(breakPoints[-1]))
#
#   output <- list(info = list(parameters = list(validFish38 = validFish38,
#                                                validBlue38 = validBlue38,
#                                                upLimitFluid120 = upLimitFluid120),
#                              n_echograms = length(breakPoints) - 1),
#                  data = data)
#
#   return(output)
# }

.getEchoData <- function(fileMode, directoryMode,
                         validFish38, validBlue38, upLimitFluid120,
                         pinInterval, date.format){

  if(is.null(fileMode) & is.null(directoryMode)){
    stop("At least whether fileMode or directoryMode must not be NULL.")
  }

  if(!is.null(fileMode)){

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

    # Clear main data (Fluid-like) using Fish and Blue noise data
    allData <- fluid120_matrix*(is.na(blue38_matrix) & is.na(fish38_matrix))
    allData[allData == 0] <- NaN
    allData <- t(allData)

  }else{

    directory <- directoryMode$directory

    # Define ttext pattern of databases
    pattern_Fish38  <- directoryMode$fish38_pattern
    # pattern_Fluid38 <- "_2Freq_Fluid38.mat"
    pattern_Blue38  <- directoryMode$blue38_pattern

    # pattern_Fish120   <- "_2Freq_Fish120.mat"
    pattern_Fluid120  <- directoryMode$fluid120_pattern
    # pattern_Blue120   <- "_2Freq_Blue120.mat"

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
      tempData <- tempData_Fluid*(is.na(tempData_Blue) & is.na(tempData_Fish))
      tempData[tempData == 0] <- NaN

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

.echogramPlot <- function(echogram, colEchogram, ...){

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

  xlim <- range(pretty_dates(xAxis, nIntervals))
  ylim <- range(pretty(yAxis, n = nIntervals))

  par(mar = c(3, 4, 2, 3), xaxs = "i", yaxs = "i")

  image(x = ext_xAxis, y = yAxis, z = newEchogram,
        xlim = .an(xlim), ylim = ylim, axes = FALSE, ylab = "Depth (m)",
        useRaster = FALSE, col = colEchogram, ...)

  axis(2, at = pretty(yAxis), labels = rev(abs(pretty(yAxis))), las = 2)
  axis(1, at = .an(pretty_dates(xlim, nIntervals)),
       labels = as.Date(pretty_dates(xlim, nIntervals)))
  axis(1, at = .an(pretty_dates(xlim, nIntervals)),
       labels = strftime(pretty_dates(xlim, nIntervals), format="%H:%M:%S"), line = 1, tick = FALSE)

  box()

  return(invisible())
}

.lineOxyrangePlot <- function(oxyrange, ...){
  xAxis <- as.POSIXct(rownames(oxyrange))

  # Add lower and upper limits of oxycline
  lines(.an(xAxis), oxyrange[,1], ...)
  lines(.an(xAxis), oxyrange[,2], ...)

  return(invisible())
}
