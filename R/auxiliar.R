.checkFilterSettings <- function(filterSettings){
  defaultFilterSettings <- get("defaultFilterSettings")

  if(is.null(filterSettings) | (is.vector(filterSettings) && is.character(filterSettings) && length(filterSettings) == 1)){
    message("Message: \nNo filter-setting object or file detected. OXim will use default filter configuration.")
    output <- defaultFilterSettings[tolower(defaultFilterSettings$name) == "default",]
  }else if(is.data.frame(filterSettings) &&
           all(is.element(c("type", "radius", "times", "tolerance"), tolower(colnames(filterSettings))))){
    output <- filterSettings
  }else{
    stop("Incorrect value for 'filterSettings'. Please define a set of filters using 'createFilterSetting' function.")
  }

  # Check variables of fileter settings object
  # Chaeck name
  if(!all(is.element(sort(unique(output$type)), c("definerFilter", "noiselessFilter"))))
    stop("Problem with 'filterSettings'. There is, at least, one wrong value on 'type' column.")

  # Check radius
  if(any(abs(as.integer(output$radius) - output$radius) > 1e-8) | any(.isOdd(output$radius)) | any(output$radius < 3))
    stop("Problem with 'filterSettings'. There is, at least, one wrong value on 'radius' column.")

  # Check times
  if(any(abs(as.integer(output$times) - output$times) > 1e-8) | any(output$times < 1))
    stop("Problem with 'filterSettings'. There is, at least, one wrong value on 'times' column.")

  # Check tolerance
  if(!is.numeric(output$tolerance) | any(output$tolerance <= 0 | output$tolerance >= 1, na.rm = TRUE))
    stop("Problem with 'filterSettings'. There is, at least, one wrong value on 'tolerance' column.")

  return(output)
}

.isOdd <- function(x){
  return(ifelse(x %% 2 != 0, FALSE, TRUE))
}

# Moving average function
movingAverage <- function(x, n = 3, circular = TRUE, ...)
{
  output <- filter(x, rep(1/n, n), circular = circular, ...)

  return(.an(output))
}

# Function for smoothing numeric vector
smoothVector <- function(x, y = NULL, ...){

  if(is.null(y)){
    y <- x
    x <- seq_along(x)
  }

  x <- .an(x)
  y <- .an(y)

  index <- !is.na(y)
  smoothFUN <- smooth.spline(x = x[index], y = y[index], ...)

  output <- predict(smoothFUN, x)

  return(output$y)
}

.getCoordsAxes <- function(coord, what){

  if(tolower(what) == "lon"){
    if(coord < 0){
      sufix <- "\u00b0 W"
    }else if(coord > 0){
      sufix <- "\u00b0 E"
    }else{
      sufix <- "\u00b0"
    }
  }else if(tolower(what) == "lat"){
    if(coord < 0){
      sufix <- "\u00b0 S"
    }else if(coord > 0){
      sufix <- "\u00b0 N"
    }else{
      sufix <- "\u00b0"
    }
  }else{
    stop("Incorrect value for 'what' parameter.")
  }

  output <- paste0(round(abs(coord), 3), sufix)

  return(output)
}

.interpIDW <- function(myData, XYZnames = c("x", "y", "z"), myGrid = NULL, ...){

  myData <- data.frame(x = myData[,XYZnames[1]],
                       y = myData[,XYZnames[2]],
                       z = myData[,XYZnames[3]],
                       stringsAsFactors = FALSE)

  coordinates(myData) <- ~ x + y

  if(is.null(myGrid)){

    xyRange <- apply(myData@coords, 2, range)

    myGrid <- expand.grid(x = seq(floor(xyRange[1, 1]), ceiling(xyRange[2, 1]), 0.01),
                          y = seq(floor(xyRange[1, 2]), ceiling(xyRange[2, 2]), 0.01),
                          stringsAsFactors = FALSE)
  }

  coordinates(myGrid) <- ~ x + y
  gridded(myGrid) <- TRUE

  myIDW <- idw(formula = z ~ 1, locations = myData, newdata = myGrid, ...)

  output <- list(myIDW = myIDW,
                 myGrid = myGrid)

  return(output)
}

# Abbreviated functions
.ac <- as.character
.an <- as.numeric
.anc <- function(...) as.numeric(as.character(...))

# Welcome message
.onAttach <- function(...) {

  packageStartupMessage("
'oXim': Tools for read Echopen outputs and get oxycline limits from echogram data.
https://cran.r-project.org/web/packages/oXim/index.html\n
'Echopen': toolbox for the multifrequency analysis of fisheries acoustics data.
http://www.france-nord.ird.fr/les-ressources/outils-informatiques\n")

  return(invisible())
}

# Get limits from cumulative sum anc
getLimits <- function(x, probs){

  x[is.na(x)] <- 0
  x <- cumsum(x)/sum(x)
  x <- cut(x = x, breaks = c(-Inf, probs, Inf), labels = c(-1, rep(0, length(probs) - 1), 1))
  x <- !abs(.anc(x))

  return(x)
}

getStartFinish <- function(x, values){

  x <- which(x)

  if(length(x) < 1){
    output <- c(1, 1)
  }else{
    output <- values[range(x)]
  }

  return(output)
}


getWeightedMatrix <- function(radius){
  constant1 <- ceiling(radius/2)
  weightedMatrix <- matrix(data = 1, nrow = radius, ncol = radius)

  for(i in 1:constant1){
    index <- seq(i, radius - (i - 1))
    weightedMatrix[i, index] <- 2
  }

  # weightedMatrix[,constant1] <- 2

  return(weightedMatrix)
}
