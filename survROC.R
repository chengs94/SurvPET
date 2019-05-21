# adapted from survivalROC

survROC <- function (Stime, status, marker, entry = NULL, predict.time,
          cut.values = NULL, method = "NNE", lambda = NULL, span = NULL,
          window = "symmetric")
{
  times = Stime
  x <- marker
  if (is.null(entry))
    entry <- rep(0, length(times))
  bad <- is.na(times) | is.na(status) | is.na(x) | is.na(entry)
  entry <- entry[!bad]
  times <- times[!bad]
  status <- status[!bad]
  x <- x[!bad]
  if (sum(bad) > 0)
    cat(paste("\n", sum(bad), "records with missing values dropped. \n"))
  if (is.null(cut.values))
    cut.values <- unique(x)
  cut.values <- cut.values[order(cut.values)]
  ncuts <- length(cut.values)
  ooo <- order(times)
  times <- times[ooo]
  status <- status[ooo]
  x <- x[ooo]
  s0 <- 1
  unique.t0 <- unique(times)
  unique.t0 <- unique.t0[order(unique.t0)]
  n.times <- sum(unique.t0 <= predict.time)
  if (method == "NNE") {
    if (is.null(lambda) & is.null(span)) {
      cat("method = NNE requires either lambda or span! \n")
      stop(0)
    }
    x.unique <- unique(x)
    x.unique <- x.unique[order(x.unique)]
    S.t.x <- rep(0, length(x.unique))
    t.evaluate <- unique(times[status == 1])
    t.evaluate <- t.evaluate[order(t.evaluate)]
    t.evaluate <- t.evaluate[t.evaluate <= predict.time]
    for (j in 1:length(x.unique)) {
      if (!is.null(span)) {
        if (window == "symmetric") {
          ddd <- (x - x.unique[j])
          n <- length(x)
          ddd <- ddd[order(ddd)]
          index0 <- sum(ddd < 0) + 1
          index1 <- index0 + trunc(n * span + 0.5)
          if (index1 > n)
            index1 <- n
          lambda <- ddd[index1]
          wt <- as.integer(((x - x.unique[j]) <= lambda) &
                             ((x - x.unique[j]) >= 0))
          index0 <- sum(ddd <= 0)
          index2 <- index0 - trunc(n * span/2)
          if (index2 < 1)
            index2 <- 1
          lambda <- abs(ddd[index1])
          set.index <- ((x - x.unique[j]) >= -lambda) &
            ((x - x.unique[j]) <= 0)
          wt[set.index] <- 1
        }
      }
      else {
        wt <- exp(-(x - x.unique[j])^2/lambda^2)
      }
      s0 <- 1
      for (k in 1:length(t.evaluate)) {
        n <- sum(wt * (entry <= t.evaluate[k]) & (times >=
                                                    t.evaluate[k]))
        d <- sum(wt * (entry <= t.evaluate[k]) & (times ==
                                                    t.evaluate[k]) * (status == 1))
        if (n > 0)
          s0 <- s0 * (1 - d/n)
      }
      S.t.x[j] <- s0
    }
    S.all.x <- S.t.x[match(x, x.unique)]
    n <- length(times)
    Sx <- sum(S.all.x[x > cut.values])/sum(x > cut.values)
  }
  return(Sx)
}
