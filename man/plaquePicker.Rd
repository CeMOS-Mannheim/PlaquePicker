% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plaquePicker.R
\name{plaquePicker}
\alias{plaquePicker}
\title{Separate stack of ion images into separated single-objects.}
\usage{
plaquePicker(
  ionImages,
  coord,
  method = c("tpoint", "geometric", "peak"),
  fixedThreshold = NULL,
  binMatrix = NULL,
  addIonImages = NULL,
  ...
)
}
\arguments{
\item{ionImages}{array of ion images as returned by \code{MALDIquant::msiSlices()}.}

\item{coord}{array of coordinates as returned by \code{MALDIquant::coordinates()}.}

\item{method}{character, method to use for binarization.
If set to "peak" the threshold is set to 0 so each signal will be counted as valid.
Make sure you set the tolerance for the ion images right
and used actual peak data to generate the ion images.}

\item{fixedThreshold}{numeric, if not NULL (default) overwrites \code{method} and uses a fixed value as a threshold.}

\item{addIonImages}{array of ion images as returned by \code{MALDIquant::msiSlices()}.
In contrast to \code{ionImages} these additional ion images will not be used for segmentation.
The data from these ion images will still be included for the summary and successive data analysis.}

\item{binMatirx}{binary matrix, if provided it overwrites \code{method} and is used to subset \code{ionImages}.
This can be helpful if you want to segment one modality (e.g. lipid data) by another (e.g. peptide data).}
}
\value{
A list, on the top level the list contains one entry per provided ion image with
    \code{conComp}      matrix, with same dimensions as the provided ion image and IDs assigned to pixels belonging to the same separate single-object.
    \code{binMat}       matrix, of ion image binarized by \code{tpoint}-function.
    \code{threshold}    numeric, the threshold estimated by \code{tpoint}-function.
    \code{spectraIdx}   list of integer vectors containing the indecies of the spectra in the list of \code{MALDIquant::massObject}.
Also, there is another entry called \code{unified}: The result of all binarized ion images combined.
    In addition to \code{conComp}, \code{binMat}, and \code{spectraIdx} that contain the equivalent matrices of the unified data as described above this entry contains the following:
    \code{intensities}  list of lists containing the intensities of the corresponding pixels for the provided ion images
}
\description{
Separate stack of ion images into separated single-objects.
}
\examples{
pp <- plaquePicker(NLGF67w_mouse1_rep1, coord = NLGF67w_mouse1_rep1_coord)
}
