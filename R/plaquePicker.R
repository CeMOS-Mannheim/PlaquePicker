#' Separate stack of ion images into seperated single-objects.
#'
#' @param ionImages array of ion images as returned by \code{MALDIquant::msiSlices()}.
#' @param coord     array of coordinates as returned by \code{MALDIquant::coordinates()}.
#'
#' @return
#' A list, on the top level the list contains one entry per provided ion image with
#'     \code{conComp}      matrix, with same dimensions as the provided ion image and IDs assigned to pixels belonging to the same seperte single-object.
#'     \code{binMat}       matrix, of ion image binarized by \code{tpoint}-function.
#'     \code{threshold}    numeric, the threshold estimated by \code{tpoint}-function.
#'     \code{spectraIdx}   list of integer vectors containing the indecies of the spectra in the list of \code{MALDIquant::massObject}.
#' Also, there is another entry called \code{unified}: The result of all binarized ion images combined.
#'     In addition to \code{conComp}, \code{binMat}, and \code{spectraIdx} that contain the equivalent matrices of the unified data as described above this entry contains the following:
#'     \code{intensities}  list of lists containing the intensities of the coresponding pixels for the provided ion images
#'
#' @export
#'
#' @examples
#' pp <- plaquePicker(NLGF67w_mouse1_rep1,
#'                          coord = NLGF67w_mouse1_rep1_coord)
plaquePicker <- function(ionImages, coord, method = c("tpoint", "geometric", "peak"), ...) {
  method <- match.arg(method)
  # helper function to extract spectra indices
  get_specIdx <- function(comp, clumpPoints) {
    spectraIdx <- vector("list", max(comp, na.rm=TRUE))
    for(clump in 1:max(comp, na.rm=TRUE) ) {
      for(pixel in 1:dim(clumpPoints[[clump]])[1]) {
        # assing MALDIquant spectra index to clump ID
        spectraIdx[[clump]][[pixel]] <- which(
          coord[,1] == as.numeric(rownames(comp)[clumpPoints[[clump]][pixel,1]])
          &
            coord[,2] == as.numeric(colnames(comp)[clumpPoints[[clump]][pixel,2]]))

      }
    }
    return(spectraIdx)
  }

  mzValues <- attr(ionImages, "center")
  resultList <- vector("list", length = length(mzValues)+1)

  # empty matrix for storage of unified binary image
  uni <- matrix(data = 0, nrow = dim(ionImages[,,1])[1],
                ncol = dim(ionImages[,,1])[2])
  # name cols and rows of uni
  # according to the original coordinate system in MALDIQuant object
  rownames(uni) <- sort(unique(coord[,1])) # x
  colnames(uni) <- sort(unique(coord[,2])) # y
  for(i in 1:length(mzValues)) {
    cat("processing", mzValues[i], i , "/", length(mzValues) ,"\n")
    # get threshold and perform binarization
    ints <- as.vector(ionImages[,,i])

    switch(method,
          "tpoint" = {
            threshold <- tpoint(ints, ...)
            cat("\n threshold of mz", mzValues[i],"based on t-point thresholding = ", threshold, "\n")
            bin <- ifelse(test = threshold < ints,
                          yes = 1,
                          no = ifelse(is.na(ints),
                                      yes = NA,
                                      no = 0))
          },
          "geometric" = {
            threshold <- geometricThreshold(ints, ...)
            cat("\n threshold of mz", mzValues[i],"based on geometric thresholding = ", threshold, "\n")
            bin <- ifelse(test = threshold < ints,
                          yes = 1,
                          no = ifelse(is.na(ints),
                                      yes = NA,
                                      no = 0))
          },
          "peak" = {
            cat("\n no threshold needed for mz", mzValues[i],"data considered as peaks\n")
            bin <- ifelse(test = ints > 0,
                          yes = 1,
                          no = ifelse(is.na(ints),
                                      yes = NA,
                                      no = 0))
            threshold <- -Inf
          }
    )




    # rebuild the matrix
    binMat <- matrix(bin,
                     nrow = dim(ionImages[,,i])[1],
                     ncol = dim(ionImages[,,i])[2])
    # perfrom connceted component labeling
    conComp <- raster::as.matrix(raster::clump(raster::raster(binMat),
                                               direction = 8))
    # name cols and rows of conComp
    # according to the original coordinate system in MALDIQuant object
    rownames(conComp) <- sort(unique(coord[,1])) # x
    colnames(conComp) <- sort(unique(coord[,2])) # y

    resultList[[i]][["conComp"]] <- conComp
    resultList[[i]][["binMat"]] <- binMat
    resultList[[i]][["threshold"]] <- threshold

    clumpPoints <- vector("list",max(conComp, na.rm=TRUE))

    # populate the list with pixels that correspont to CC ID
    cat("extracting clump information...\n")
    for (j in 1:max(conComp, na.rm=TRUE))
    {
      clumpPoints[[j]]            = which(conComp == j, arr.ind = TRUE)
    }

    spectraIdx <- get_specIdx(comp = conComp,
                              clumpPoints = clumpPoints)
    resultList[[i]][["spectraIdx"]] <- spectraIdx

    # add binMat to unified binary image
    uni <- uni + binMat
    uni[uni > 0] <- 1
  }

  resultList[[length(resultList)]][["uniBin"]] <- uni
  uniComp <- raster::as.matrix(raster::clump(raster::raster(uni),
                                             direction = 8))

  # name cols and rows of uniComp
  rownames(uniComp) <- rownames(uni)
  colnames(uniComp) <- colnames(uni)
  resultList[[length(resultList)]][["uniComp"]] <- uniComp

  clumpPoints <- vector("list",max(uniComp, na.rm=TRUE))

  # populate the list with pixels that correspont to CC ID
  cat("extracting unified clump information...\n")
  for (i in 1:max(uniComp, na.rm=TRUE))
  {
    clumpPoints[[i]]            = which(uniComp == i, arr.ind = TRUE)
  }

  spectraIdx <- get_specIdx(comp = uniComp,
                            clumpPoints = clumpPoints)
  resultList[[length(resultList)]][["spectraIdx"]] <- spectraIdx
  names(resultList) <- c(mzValues, "unified")
  intensities <- vector("list", max(uniComp, na.rm=TRUE))
  for(i in 1:length(intensities)) {
    extractedInt <- vector("list", length(mzValues))
    names(extractedInt) <- mzValues
    for(j in 1:length(mzValues)) {

      extractedInt[[j]] <-ionImages[,,j][which(0<(uniComp==i))]
    }
    intensities[[i]] <- extractedInt
  }
  resultList[[length(resultList)]][["intensities"]] <- intensities
  return(resultList)
}
