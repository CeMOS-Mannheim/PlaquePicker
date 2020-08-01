#' Quantile correction
#' @description
#' Function to perform quantile correction for better visual representation of ion images.
#'
#' @param im matrix containing intensities of ion image.
#' @param q  quantile that is set to 100% intensity
#'
#' @return
#' matrix of quantile corrected intensities.
#' @export
#'
#' @examples
#' im_cor <- qcor(NLGF67w_mouse1_rep1[,,1], 0.9995)

qcor <- function(im, q = 0.999) {
  im_cor <- im/quantile(im, q, na.rm = T)*100
  im_cor[which(im_cor>100)] <- 100
  return(im_cor)
}
