#' get mz values in plaquePicker result list
#'
#' @param pp plaquePicker result list
#'
#' @return
#' vector of mz values as character.
#' @export
get_mzValues <- function(pp) {
  nm <- names(pp)
  return(nm[-length(nm)])
}
