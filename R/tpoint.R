#' T-Point thresholding
#'
#' @description
#' implemation described in Coudray et. al. 2010,
#' see https://www.sciencedirect.com/science/article/pii/S0167865509003675?via%3Dihub#fd6
#'
#' @param ints      numeric vector of intensities.
#' @param breaks    numeric, number of bins in histogram.
#' @param diagnosis logical, return set of internal values additional to threshold
#' @param plot      logical, show plot for diagnosis.
#'
#' @return
#' numeric threshold.
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_col geom_line geom_vline theme_bw geom_text labs theme
#' @importFrom graphics hist
#'
#' @examples
#' thresh <- tpoint(as.vector(NLGF67w_mouse1_rep1[,,1]))
tpoint <- function(ints, breaks = 500, diagnosis = FALSE, plot = FALSE) {

  # remove zero values
  ints <- as.vector(ints[which(!ints==0)])

  hist <- hist(ints, breaks = breaks, plot = FALSE)

  bin_val_full <- hist$mids
  counts_full <- hist$counts
  M = which.max(counts_full)
  L = length(counts_full)

  n   = 0
  Sh  = 0
  Shh = 0
  Sg  = 0
  Sgg = 0
  Sgh = 0
  e1 = vector("numeric", length = L)

  for(k in M:L) {
    n = n + 1
    Sh = sum(counts_full[M:k])
    Shh = sum(counts_full[M:k]^2)
    Sg = sum(bin_val_full[M:k])
    Sgg = sum(bin_val_full[M:k]^2)
    Sgh = sum(counts_full[M:k] * bin_val_full[M:k])
    e1[k] = Shh - Sh/n - (n * Sgh - Sg * Sh)^2/(n * (n * Sgg - Sg^2))
  }

  n   = 0
  Sh  = 0
  Shh = 0
  Sg  = 0
  Sgg = 0
  Sgh = 0
  e2 = vector("numeric", length = L)
  e2[L] = 0

  for(k in L:(M+1)) {
    n = n + 1
    Sh = sum(counts_full[k:L])
    Shh = sum(counts_full[k:L]^2)
    Sg = sum(bin_val_full[k:L])
    Sgg = sum(bin_val_full[k:L]^2)
    Sgh = sum(counts_full[k:L] * bin_val_full[k:L])
    e2[(k-1)] = Shh - Sh/n - (n * Sgh - Sg * Sh)^2/(n * (n * Sgg - Sg^2))
  }
  tpoint = which.min(e1[M:L] + e2[M:L])

  df_list <- list("hist" = hist,
                  "eValues" = tibble::tibble(
                    breaks = hist$breaks[-1],
                    counts = hist$counts,
                    e1 = e1,
                    e2 = e2,
                    sumE =e1 + e2),
                  "res"= tibble::tibble(M = M,
                                        L = L,
                                        tpoint_bin = M+tpoint,
                                        threshold = bin_val_full[M+tpoint]))

  if(plot) {
    counts_df <- tibble(breaks = df_list$hist$breaks[-1], counts = df_list$hist$counts)
    df <- tibble(breaks = df_list$hist$breaks[-355],
                 counts = df_list$hist$counts,
                 "Error e1 (near peak)" = df_list$eValues$e1,
                 "Error e2 (near tail)" = df_list$eValues$e2,
                 "Error sum" = df_list$eValues$sumE)

    df$"Error e1 (near peak)"[which(is.na(df$"Error e1 (near peak)"))] <- 0

    df <- df %>%
      gather(measure, value, -breaks, -counts) %>%
      group_by(measure) %>%
      mutate(value = value/max(value, na.rm = TRUE))


    ggplot(counts_df, aes(x = breaks, y = counts/max(counts))) +
      geom_col(col = "grey75") +
      geom_line(data = df, aes(y = value, col = measure), size = 1) +
      geom_vline(aes(xintercept =df_list$hist$breaks[df_list$res$tpoint_bin]), linetype = 2, size = 1.2) +
      theme_bw() +
      geom_text(data = df_list$res, aes(x = threshold + threshold * 0.25, y = 0.5, label =  round(threshold,4))) +
      theme(text = element_text(size = 12),
            legend.text = element_text(face = "bold"),
            legend.position = "bottom") +
      labs(col = "",
           x = "Intensity",
           y = "rel. Counts") -> p
    print(p)
  }

  if(diagnosis) {

    return(df_list)
  }



  return(bin_val_full[M+tpoint])
}

