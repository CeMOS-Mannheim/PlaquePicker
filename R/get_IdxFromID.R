#' Get indicies of spectra per plaque ID
#'
#' @param pp       plaque picker result list. See \code{plaquePicker}
#' @param ID       vector of ID, set to NA to return all Indicies associated with any clump.
#' @param mzValue  mzValue to ID is refering to. By default the unified IDs are used but this can also be overwritten by entering a specific mz.
#'
#' @return
#' vector of indicies refering to the spectra in the MALDIquant spectra list used to construct the ion images used.
#' @export

get_IdxFromID <- function(pp, ID, mzValue = "unified") {
  mzValue <- as.character(mzValue)
  if(length(mzValue) > 1) {
    stop("mzValue needs to be length == 1\n")
  }
  if(!mzValue %in% c(get_mzValues(pp), "unified")) {
    stop("mzValue not contained in pp\n")
  }
  if(any(is.na(ID))) {
    cat("ID set to NA. Returning all indicies associated with any clump.\n")
    ID <- 1:length(pp[[mzValue]][["spectraIdx"]])
  }
  idx <- lapply(ID, function(x) {
    pp[[mzValue]][["spectraIdx"]][[x]]
  })
  return(unlist(idx))
}
