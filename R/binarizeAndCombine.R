#' Binarize images and combine to a single binary image
#'
#' @param im          matrix or array representing a image or a stack of images
#' @param method      character, can be "tpoint", "geometric", "peak"
#' @param binMatrix   matrix, a binary matrix. Overwrites method and uses this matrix to binarize the data.
#'
#' @return
#' matrix with binary image
#' @export

binarizeAndCombine <- function(im, method = c("tpoint", "geometric", "peak"), binMatrix = NULL, ...) {
  if(!is.null(binMatrix)) {
    # check if dimof ionImages and binMatrix match. Otherwise stop.
    if(!dim(im)[1] == dim(binMatrix)[1] & !dim(im)[2] == dim(im)[2]) {
      stop("Dimensions of binMatrix has to be the same as ionImages! \n
           !dim(ionImages[,,1]) == dim(binMatrix)\n")
    }
    method <- "binMat"
  } else {
    method <- match.arg(method)
  }

  if(is.matrix(im)) {
    z <- 1
  } else {
    z <- dim(im)[3]
  }

  uni <- matrix(data = 0, nrow = dim(im)[1],
                ncol = dim(im)[2])

  for(i in 1:z) {
    if(is.matrix(im)) {
      ints <- as.vector(im)
    } else {
      ints <- as.vector(im[,,i])
    }

    switch(method,
           "tpoint" = {
             threshold <- tpoint(ints, ...)
             bin <- ifelse(test = threshold < ints,
                           yes = 1,
                           no = ifelse(is.na(ints),
                                       yes = NA,
                                       no = 0))
           },
           "geometric" = {
             threshold <- geometricThreshold(ints, ...)
             bin <- ifelse(test = threshold < ints,
                           yes = 1,
                           no = ifelse(is.na(ints),
                                       yes = NA,
                                       no = 0))
           },
           "peak" = {
             bin <- ifelse(test = ints > 0,
                           yes = 1,
                           no = ifelse(is.na(ints),
                                       yes = NA,
                                       no = 0))
             threshold <- -Inf
           },
           "binMat" = {
             bin <- binMatrix
             threshold <- -Inf
           }
    )

    # rebuild the matrix
    binMat <- matrix(bin,
                     nrow = dim(im)[1],
                     ncol = dim(im)[2])
    uni <- uni + binMat
    uni[uni > 0] <- 1
  }
  return(uni)
  }
