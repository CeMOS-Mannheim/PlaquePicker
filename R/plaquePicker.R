#' Internal helper function to extract spectra indices
#'
#' @param comp        matrix, result of \code{raster::clump}. Matrix with indicies for each unique connected component.
#' @param coord       array of coordinates as returned by \code{MALDIquant::coordinates()}.
#'
#' @return
#' list of spectra indicies.
get_specIdx <- function(comp, coord) {
  maxComp_idx <- max(comp, na.rm=TRUE)

  clumpPoints <- vector("list", maxComp_idx)
  # populate the list with pixels that correspont to CC ID
  for (j in 1:maxComp_idx)
  {
    clumpPoints[[j]]            = which(comp == j, arr.ind = TRUE)
  }

  coord_x <- coord[,1]
  coord_y <- coord[,2]

  x <- rownames(comp)
  y <- colnames(comp)

  spectraIdx <- vector("list", maxComp_idx)
  for(clump in 1:maxComp_idx) {
    numPx <- dim(clumpPoints[[clump]])[1]
    for(pixel in 1:numPx) {
      # assign MALDIquant spectra index to clump ID
      spectraIdx[[clump]][[pixel]] <- which(
        coord_x == as.numeric(x[clumpPoints[[clump]][pixel,1]])
        &
          coord_y == as.numeric(y[clumpPoints[[clump]][pixel,2]]))
    }
  }
  return(spectraIdx)
}

#' Get intensity values from ion images
#'
#' @param comp matrix, locations of connected components
#' @param ii   array, ion images
#'
#' @return
#' list of intensities.
get_intensities <-   function(comp, ii) {
  mzValues <- attr(ii, "center")
  intensities <- vector("list", max(comp, na.rm=TRUE))
  for(i in 1:length(intensities)) {
    extractedInt <- vector("list", length(mzValues))
    names(extractedInt) <- mzValues
    for(j in 1:length(mzValues)) {

      extractedInt[[j]] <-ii[,,j][which(0<(comp==i))]
    }
    intensities[[i]] <- extractedInt
  }
  return(intensities)
}



#' Separate stack of ion images into separated single-objects.
#'
#' @param ionImages array of ion images as returned by \code{MALDIquant::msiSlices()}.
#' @param coord     array of coordinates as returned by \code{MALDIquant::coordinates()}.
#' @param method    character, method to use for binarization.
#'                  If set to "peak" the threshold is set to 0 so each signal will be counted as valid.
#'                  Make sure you set the tolerance for the ion images right
#'                  and used actual peak data to generate the ion images.
#' @param binMatirx binary matrix, if provided it overwrites \code{method} and is used to subset \code{ionImages}.
#'                  This can be helpful if you want to segment one modality (e.g. lipid data) by another (e.g. peptide data).
#'
#' @return
#' A list, on the top level the list contains one entry per provided ion image with
#'     \code{conComp}      matrix, with same dimensions as the provided ion image and IDs assigned to pixels belonging to the same separate single-object.
#'     \code{binMat}       matrix, of ion image binarized by \code{tpoint}-function.
#'     \code{threshold}    numeric, the threshold estimated by \code{tpoint}-function.
#'     \code{spectraIdx}   list of integer vectors containing the indecies of the spectra in the list of \code{MALDIquant::massObject}.
#' Also, there is another entry called \code{unified}: The result of all binarized ion images combined.
#'     In addition to \code{conComp}, \code{binMat}, and \code{spectraIdx} that contain the equivalent matrices of the unified data as described above this entry contains the following:
#'     \code{intensities}  list of lists containing the intensities of the corresponding pixels for the provided ion images
#'
#' @export
#'
#' @examples
#' pp <- plaquePicker(NLGF67w_mouse1_rep1, coord = NLGF67w_mouse1_rep1_coord)
plaquePicker <- function(ionImages, coord, method = c("tpoint", "geometric", "peak"), binMatrix = NULL, ...) {
  if(!is.null(binMatrix)) {
    if(!dim(ionImages[,,1]) == dim(binMatrix)) {
      stop("Dimensions of binMatrix has to be the same as ionImages! \n
           !dim(ionImages[,,1]) == dim(binMatrix)\n")
    }
    method <- "binMat"
  } else {
    method <- match.arg(method)
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
           },
           "binMat" = {
             cat("\n no threshold needed for mz", mzValues[i],"binMatrix used for segmentation\n")
             bin <- binMatrix
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

    cat("extracting clump information...\n")
    resultList[[i]][["spectraIdx"]] <- get_specIdx(comp = conComp, coord = coord)

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

  cat("extracting unified clump information...\n")
  resultList[[length(resultList)]][["spectraIdx"]] <- get_specIdx(comp = uniComp, coord = coord)
  names(resultList) <- c(mzValues, "unified")


  resultList[[length(resultList)]][["intensities"]] <- get_intensities(comp = uniComp, ii = ionImages)
  return(resultList)
}
