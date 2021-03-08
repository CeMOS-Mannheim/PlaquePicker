#' Combine two stacks ion images into one stack of ion images
#'
#' @description
#' This function can be used to combine different stacks of ion images (e.g. from different polarities).
#' For this to work properly the images have to be aligned.
#' Small differences (e.g. one additional row) can be compensated with this function.
#' Be sure you know what you are doing and check if the quality of the result is sufficent for your means!
#'
#' @param ii1 Array, ion images as returend by MSISclices from MALDIquant
#' @param ii2 Array, ion images as returend by MSISclices from MALDIquant
#' @param add Logical, if the dimensions of ii1 and ii2 dont match, add rows and cols to compensate
#'
#' @return
#' A stack of combined ion images
#' @export

combineIonImages <- function(ii1, ii2, add = TRUE) {
  d1 <- dim(ii1)
  d2 <- dim(ii2)
  if(!length(d1) == 3 | !length(d2) == 3) {
    # check if both dims are 3-dimenional
    stop("ii1 or ii2 does not have 3 dimensions!\n")
  }
  mz1 <- attr(ii1, "center")
  mz2 <- attr(ii2, "center")

  if(!all(d1[1:2] == d2[1:2])) {
    # check if x and y dim match
    if(!add) {
      stop("dimensions of ii1 and ii2 do not match!\n")
    }
    warning("dimensions of ii1 and ii2 do not match! Adding rows and cols to compensate.\n")
    dX <- c(d1[1], d2[1])
    maxX <- dX[which.max(dX)]
    dY <- c(d1[2], d2[2])
    maxY <- dY[which.max(dY)]


    ii1New <- array(dim = c(x = maxX, y = maxY, z = d1[3]))
    ii2New <- array(dim = c(x = maxX, y = maxY, z = d2[3]))

    # correct ii1
    if(!d1[1] == maxX | !d1[2] == maxY) {
      for(z in 1:dim(ii1New)[3]) {
        for(x in 1:dim(ii1New)[1]) {
          diffY <- abs(d1[2]-maxY)
          if(x <= dim(ii1)[1]) {
            newVec <- c(as.vector(ii1[x,,z]),
                        rep(NA, diffY))
          } else {
            newVec <- rep(NA, maxY)
          }
          ii1New[x,,z] <- newVec
        }
      }
      ii1 <- ii1New
    }
    if(!d2[1] == maxX | !d2[2] == maxY) {
      for(z in 1:dim(ii2New)[3]) {
        for(x in 1:dim(ii2New)[1]) {
          diffY <- abs(d2[2]-maxY)
          if(x <= dim(ii2)[1]) {
            newVec <- c(as.vector(ii2[x,,z]),
                        rep(NA, diffY))
          } else {
            newVec <- rep(NA, maxY)
          }
          ii2New[x,,z] <- newVec
        }
      }
      ii2 <- ii2New
    }
    # update dimensions for new array
    d1 <- dim(ii1)
    d2 <- dim(ii2)
  }



  iiAll <- array(data = c(ii1, ii2),
                 dim = c(d1[1:2], d1[3]+d2[3]))
  attr(iiAll, "center") <- c(mz1, mz2)
  return(iiAll)
}
