#' Local FSR to svalue
#'
#' Convert local FSR to svalue by taking the cumulative mean of
#' increasing local FSRs.
#' This is used within \code{apeglm} to generate the svalue table,
#' but provided as convenience function.
#'
#' @param lfsr the local false sign rates
#'
#' @return s-values
#'
#' @references
#'
#' False sign rate and s-value:
#'
#' Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2. doi: 10.1093/biostatistics/kxw041
#' 
#' @export
#' 
#' @examples 
#' 
#' # Example 1 -- simulated local FSR data
#' local.fsr <- runif(1000)
#' sval <- svalue(local.fsr)
#'
#' # Example 2 -- first runs example from 'apeglm'
#' example("apeglm")
#' local.fsr <- res$fsr[1, ]
#' sval <- svalue(local.fsr)
#' 
svalue <- function(lfsr) {
  lfsr.sorted <- sort(lfsr, na.last = TRUE)
  cumsum.idx <- seq_along(lfsr)
  cumsum.lfsr <- cumsum(lfsr.sorted)
  (cumsum.lfsr/cumsum.idx)[rank(lfsr, ties.method = "first", na.last = TRUE)]
}
