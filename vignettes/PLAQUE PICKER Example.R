## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi=600,fig.width=8)


## -----------------------------------------------------------------------------
install.packages("MALDIquant", dependencies = TRUE)
install.packages("MALDIquantForeign", dependencies = TRUE)

# install devtools to be able to install packages from github
install.packages("devtools", dependencies = TRUE)
devtools::install_github("CeMOS-Mannheim/plaquePicker")

# The following dependencies are not needed to use the core functions of the package
# but will be used in this vignette
install.packages("tidyverse", dependencies = TRUE)
install.packages("viridis", dependencies = TRUE)

# for the plotting of Venn-Diagrams we will be using the package Vennerable
# which is not in the CRAN repository and has some dependencies to the
# BioConductor repository.
if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("RBGL")
BiocManager::install("graph")
devtools::install_github("js229/Vennerable")


## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(PlaquePicker)
  library(MALDIquant)        # general MS functions
  library(MALDIquantForeign) # for import of imzML data
  library(tidyverse)         # general data science tools
  library(viridis)           # pretty colors for ion images

  # for the plotting of Venn-Diagrams we will be using the package Vennerable
  # which is not in the CRAN repository and has some dependencies to the
  # BioConductor repository. To install Vennerable and its dependencies
  # use the following commands:
  #
  # if (!requireNamespace("BiocManager", quietly = TRUE))
  #     install.packages("BiocManager")
  # BiocManager::install(version = "3.11")
  # BiocManager::install("RBGL")
  # BiocManager::install("graph")
  # devtools::install_github("js229/Vennerable")
  library(Vennerable)})


## ---- eval = TRUE, results='hide'---------------------------------------------
# unzip spectra
unzip("data-raw/NLGF67w_mouse1_rep1.zip")

# read spectra
spec <- importImzMl("NLGF_Prot_NLGF1.imzML",
                    verbose = FALSE)

# tidy up
file.remove("NLGF_Prot_NLGF1.imzML")
file.remove("NLGF_Prot_NLGF1.ibd")



## ---- warning=FALSE, eval=FALSE-----------------------------------------------
## spec <- calibrateIntensity(spec,
##                            method = "TIC")
##
## # small halfWindowSize needed as number of points
## # per spectra reduced for smaller example data sets
## spec <- smoothIntensity(spec,
##                         method = "SavitzkyGolay",
##                         halfWindowSize = 2)
##
## spec <- removeBaseline(spec,
##                        method = "TopHat")
##
##


## -----------------------------------------------------------------------------
avgSpec <- averageMassSpectra(spec)

# shorten file path (for plotting reasons only)
avgSpec@metaData$file <- basename(avgSpec@metaData$file)

plot(avgSpec, ylab = "Intensity [a.u.]")
lines(detectPeaks(avgSpec))
labelPeaks(detectPeaks(avgSpec),
           digits = 0)


## -----------------------------------------------------------------------------
ionImages <- msiSlices(spec,
                       center = c(4059.9,
                                  4159.1,
                                  4442.6),
                       tolerance = 5)


## -----------------------------------------------------------------------------
par(mfrow = c(1, 3), mar =c(0,0,0,0))
image(qcor(ionImages[,,1], 0.9995),
      col = viridis::cividis(30),
      asp = 1,
      axes = FALSE)
title("Ab1-38Arc", line = -2)
image(qcor(ionImages[,,2], 0.9995),
      col = viridis::cividis(30),
      asp = 1,
      axes = FALSE)
title("Ab1-39Arc", line = -2)
image(qcor(ionImages[,,3], 0.9995),
      col = viridis::cividis(30),
      asp = 1,
      axes = FALSE)
title("Ab1-42Arc", line = -2)


## -----------------------------------------------------------------------------
hist(ionImages[,,1],
     breaks = 300,
     xlab = "Intensity [a.u.]",
     main = "Histogram of Ab1-38Arc intensity")
thresh_Ab38 <- tpoint(ionImages[,,1], plot = TRUE)


## -----------------------------------------------------------------------------
bin <- ifelse(test = thresh_Ab38 < as.vector(ionImages[,,1]),
              yes = 1,
              no = ifelse(is.na(as.vector(ionImages[,,1])),
                          yes = NA,
                          no = 0))
# rebuild the matrix
binMat <- matrix(bin,
                 nrow = dim(ionImages[,,1])[1],
                 ncol = dim(ionImages[,,1])[2])
par(mfrow = c(1, 2), mar = c(0,0,0,0))
image(qcor(ionImages[,,1], 0.9999),
      col = viridis::cividis(30),
      asp = 1,
      axes = FALSE)
title("Ab1-38Arc intensities", line = -1)
image(binMat,
      col = c("black", "white"),
      asp = 1,
      axes = FALSE)
title("Ab1-38Arc binarized", line = -1)


## ---- warning=FALSE-----------------------------------------------------------
coord <- coordinates(spec)
pp <- plaquePicker(ionImages = ionImages,
                   coord = coord,
                   method = "tpoint")


## -----------------------------------------------------------------------------
# first extract the MALDIquant indicies associated with plaque
idx <- get_IdxFromID(pp, ID = NA) # setting ID to NA will extract all IDs instead of specific IDs

# compute average spectrum of plaque associated pixels
plaqueAvg <- averageMassSpectra(spec[idx])

# shorten file path (for plotting reasons only)
plaqueAvg@metaData$file <- basename(plaqueAvg@metaData$file)

plot(plaqueAvg,
     ylab = "Intensity [a.u.]")
lines(avgSpec,
      lty=2)
labelPeaks(detectPeaks(plaqueAvg),
           digits = 0)
lines(detectPeaks(plaqueAvg))
legend("right",
       legend=c("Plaque avg. spectrum","Overall avg. spectrum"),
       lty=1:2,
       cex=0.8)


## -----------------------------------------------------------------------------
plot_venn(pp, mzNames = c("Ab1-38Arc",
                          "Ab1-39Arc",
                          "Ab1-42Arc"),
          relative = TRUE)


## -----------------------------------------------------------------------------
df <- bind_rows(pp$unified$intensities, .id = "ID") %>%
        mutate(ID = as.numeric(ID)) %>%
        group_by(ID) %>%
        mutate(size = n() * 400)
head(df)


## -----------------------------------------------------------------------------
ggplot(df, aes(x = size)) +
  geom_histogram(bins = 40) +
  theme_bw() +
  labs(x = "Plaque area [um²]")



## -----------------------------------------------------------------------------
df_mean <- df %>%
  ungroup() %>%
  rename("Ab1_38Arc" = "4059.9",
         "Ab1_39Arc" = "4159.1",
         "Ab1_42Arc" = "4442.6") %>%
  group_by(ID) %>%
  gather(mz, int, -size, -ID) %>%
  group_by(ID, mz) %>%
  summarise(meanInt = mean(int),
            size = first(size),
            .groups = "drop_last")
df_mean %>%
  spread(mz, meanInt) %>%
  head()



## -----------------------------------------------------------------------------
ggplot(df_mean, aes(x = mz, y = meanInt)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Abeta species",
       y = "Mean plaque Int. [a.u.]")




## -----------------------------------------------------------------------------
df_size <- df_mean %>%
  spread(mz, meanInt) %>%
  group_by(ID) %>%
  mutate(sizeGroup = ifelse(size <= 400,
                            "small",
                            ifelse(size <= 2000,
                                   "medium",
                                   "big")))

df_size %>%
  group_by(sizeGroup) %>%
  summarise(n = length(unique(ID)), .groups = "drop_last") %>%
  ggplot(aes(x = sizeGroup,
             y = n,
             label = n)) +
  geom_col() +
  geom_text(y = 100) +
  theme_bw() +
  labs(x = "Size group",
       y = "Number of observations")




## -----------------------------------------------------------------------------
df_size %>%
  mutate(ratio = Ab1_42Arc/Ab1_38Arc) %>%
  ggplot(aes(x = sizeGroup, y = ratio)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Size Group",
       y = "Ratio Ab1-42Arc/Ab1-38Arc")


## -----------------------------------------------------------------------------
sessionInfo()

