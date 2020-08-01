#' point2LineDist
#'
#' @description
#' Calculates the distance between point \code{x} to the line defined by \code{p1} and \code{p2}.
#' ref: https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
#'
#' @param p1 Point 1 c(x,y)
#' @param p2 Point 2 c(x,y)
#' @param x  Point which distance should be calculated in relation to line formed by \code{p1} and \code{p2}.
#'
#' @return
#' Distance of \code{x} to line defined by \code{p1} and \code{p2}

point2LineDist      = function(p1, p2, x) {
  # suppose you have a line defined by two points P1 and P2, then the
  # distance from point x to the line is defined (in 2D) by

  d         = abs(((p2[2] - p1[2]) * x[1]) - ((p2[1] - p1[1]) * x[2]) + (p2[1] * p1[2]) - (p2[2] * p1[1]) /
                    sqrt(sum((p2 - p1) ** 2)))

  return(d)

}

#' Unimodal threshold based on geometric method
#'
#' @description
#' Calculates a threshold for unimodal distribuated values based on histogram thresholding.
#' Implementation of Rosin, P. L. (2001). "Unimodal thresholding." Pattern Recognition 34(11): 2083-2096.
#' A line is drawn between the bin with the highest frequency and the bin with the highest intensity/value (main line).
#' After for each bin another line, orthogonal to the first line, is drawn to the "top" of each bin.
#' The value of the bin with the longst distance to the main line is defined as threshold
#'
#' @param intensities Vector of intensities.
#' @param breaks            Number of bins of the histogram.
#' @param plot          Plots the internally generated plot.
#'
#' @return
#' Threshold value (see description)
#'
#'
#' @export

geometricThreshold <- function(intensities, breaks = 500, plot = FALSE) {


  # Generate histogram object
  H                   = hist(intensities, breaks = breaks, plot = FALSE)

  # find the coordinates of the max bin (peak) and the last bin (tail)
  pIdx                = which.max(H$counts)
  tIdx                = length(H$mids)

  peak                = c(H$mids[pIdx], H$counts[pIdx])
  tail                = c(H$mids[tIdx], H$counts[tIdx])

  if(plot) { # Cosmetic: plots dots on peak and tail and connects them by a red line
    hist(intensities, breaks = breaks,
         xlab = "Intensities",
         ylab = "Frequency")
    points(x = peak[1], y = peak[2], col = "red", pch = 20, cex = 2)
    points(x = tail[1], y = tail[2], col = "red", pch = 20, cex = 2)
    lines(x = rbind(peak, tail), col = "red", lty = 2)
  }

  # Generate empty vector
  allDist             = vector(mode = "numeric", length = length(H$mids) - pIdx + 1)

  # loop through each bin and compute the distance to the line connecting the peak to the tail
  for (ii in pIdx : length(H$mids))
  {

    # compute the distance
    allDist[ii]         = point2LineDist(peak, tail, c(H$mids[ii], H$counts[ii]))



  }

  # find the max distance
  maxIdx              = which.max(allDist)

  if(plot) {

    points(x = H$mids[maxIdx], y = H$counts[maxIdx], col = "green", pch = 10, cex = 5)
  }

  # Finally threshold your intensityObject based on the resulting threshold.
  # end
  return(H$mids[maxIdx])
}

