#' Generate a Venn-diagram for a plaquePicker result list
#'
#' @param pp       plaquePicker result list
#' @param mzIdx    integer vector, indicies of mz values. To get an overview of mz values try \code{get_mzValues()}. Can not be longer then 3.
#' @param mzNames  character vector, names of mz values (e.g. "ABeta1-38"). If set to NULL (default) the mz values themself will be used as names during plotting.
#' @param plot     logical, if set to FALSE the \code{Venn} object is returned.
#' @param method   character, if set to "plaque" the Venn diagram will be based on overlapping signals per plaque. If set to "pixel" then it will be pixel-wise.
#' @param relative logical, if set to false the number of plaques is returned instead of the percentage.
#'
#'
#' @export
#'
#' @importFrom stats na.omit
#' @importFrom magrittr %>%
plot_venn <- function(pp, mzIdx = 1:3, mzNames =NULL, plot = TRUE, method = c("plaque", "pixel"), relative = TRUE) {
  method = match.arg(method)
  if(!requireNamespace("Vennerable")) {
    stop("To use this function the package 'Vennerable' is needed.\n
         Install it using devtools::install_github('js229/Vennerable').\n
         Vennerable depends on 'graph' and 'RBGL' from Bioconductor!\n")
  }
  if(!is.null(mzNames)) {
    if(!length(mzNames) == length(mzIdx)) {
      stop("mzNames must be either NULL or same lenght as mzIdx\n")
    }
  }
  # extract matrix of unified plaque IDs
  unified <- pp$unified$uniComp
  # get mzValues of IonImages
  mzVal <- get_mzValues(pp)[mzIdx]
  PlaqueIDs_venn <- vector("list",
                           length = length(mzIdx))
  for(i in 1:length(mzIdx)) {
    switch(method,
           "plaque" = {
             PlaqueIDs_venn[[i]] <- pp[[mzIdx[i]]]$binMat * unified
             PlaqueIDs_venn[[i]] <- PlaqueIDs_venn[[i]] %>%
               as.vector() %>%
               unique() %>%
               na.omit() %>%
               sort() %>%
               .[-1]
           },
           "pixel" = {
             PlaqueIDs_venn[[i]] <- unlist(pp[[mzIdx[i]]]$spectraIdx)
           })
  }
  if(!is.null(mzNames)) {
    names(PlaqueIDs_venn) <- mzNames
  } else {
    names(PlaqueIDs_venn) <- mzVal
  }
  venn <- Vennerable::Venn(PlaqueIDs_venn)

  if(relative) {
    l <- length(mzIdx)+1
    venn@IndicatorWeight[,l] <-
      round(venn@IndicatorWeight[,l]/sum(venn@IndicatorWeight[,l]),3) * 100
  }
  if(!plot) {
    return(venn)
  }
  if(length(mzIdx) < 4) {
  Vennerable::plot(venn)
  } else {
    Vennerable:::plotVenn(Vennerable:::compute.S4(venn))
  }
}
