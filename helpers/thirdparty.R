deseq2_prevent_aggressive_filtering = function(numRej, theta) {
    # Copyright DEseq2 authors: Michael Love, Simon Anders, Wolfgang Huber
    # Licence: LGPL (>= 3)
    # Source file: https://github.com/mikelove/DESeq2/blob/master/R/results.R
    lo.fit <- lowess(numRej ~ theta, f=1/5)
    if (max(numRej) <= 10) {
      j <- 1
    } else { 
      residual <- if (all(numRej==0)) {
        0
      } else {
        numRej[numRej > 0] - lo.fit$y[numRej > 0]
      }
      thresh <- max(lo.fit$y) - sqrt(mean(residual^2))
      j <- if (any(numRej > thresh)) {
        which(numRej > thresh)[1]
      } else {
        1 
      }
    }
    j
}
