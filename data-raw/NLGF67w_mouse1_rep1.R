## code to prepare `NLGF67w_mouse1_rep1`
## collection of ion images for Abeta species Ab1-38Arc (4059.9), Ab1-39Arc (4159.1) and Ab1-42Arc (4442.6) with 5 Da tolerance each.

# the zip folder contains a MSI image of a full brain (coronal section) of a NLGF mouse
# at 67 weeks of age at bregma -1.5mm in the mz rage of 4k-4.5k and reduced spectral resolution (1/4)


# unzip spectra
unzip("data-raw/NLGF67w_mouse1_rep1.zip")

# read spectra
NLGF67w_mouse1_rep1_spec <- MALDIquantForeign::importImzMl("NLGF_Prot_NLGF1.imzML",
                                       verbose = FALSE)

# preprocess spectra
NLGF67w_mouse1_rep1_spec <- MALDIquant::calibrateIntensity(NLGF67w_mouse1_rep1_spec,
                                       method = "TIC")
NLGF67w_mouse1_rep1_spec <- MALDIquant::smoothIntensity(NLGF67w_mouse1_rep1_spec,
                                    method = "SavitzkyGolay",
                                    halfWindowSize = 2) # small halfWindowSize needed as number of points per spectra reduced for smaller example datasets
NLGF67w_mouse1_rep1_spec <- MALDIquant::removeBaseline(NLGF67w_mouse1_rep1_spec,
                                   method = "TopHat")
NLGF67w_mouse1_rep1_coord <- MALDIquant::coordinates(NLGF67w_mouse1_rep1_spec)
NLGF67w_mouse1_rep1 <- MALDIquant::msiSlices(NLGF67w_mouse1_rep1_spec,
                                             center = c(4059.9,
                                                        4159.1,
                                                        4442.6),
                                             tolerance = 5)

# tidy up
file.remove("NLGF_Prot_NLGF1.imzML")
file.remove("NLGF_Prot_NLGF1.idb")

# store processed data for later use as example data
usethis::use_data(NLGF67w_mouse1_rep1, overwrite = TRUE)
usethis::use_data(NLGF67w_mouse1_rep1_coord, overwrite = TRUE)
