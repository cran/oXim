#' @rdname echogramPlot
#' @export
echogramPlot.default <- function(x, colEchogram = "colPalette", ...){
  colEchogram <- get(colEchogram)
  .echogramPlot(echogram = x, colEchogram = colEchogram, ...)

  return(invisible())
}
