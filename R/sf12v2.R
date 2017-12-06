#' SF12v2 questionnaire scoring
#'
#' SF12v2 questionnaire scoring
#' @param X a \code{\link{matrix}} or \code{\link{data.frame}} of 12
#' columns, containing questionnaire items. In order from left to right:
#' gh1, pf02, pf04, rp2, rp3, re2, re3, bp2, mh3, vt2, mh4, sf2.
#' @note
#' This is an R port of SAS algorithm by Apolone and Mosconi found
#' \href{http://crc.marionegri.it/qdv/index.php?page=sf12}{here}.
#' 
#' SF-12  is a registered trademark of medical outcomes trust.
#' @examples
#' ## -------------------------
#' ## Algorithm test/validation
#' ## -------------------------
#' (scores <- sf12(sf12sample))
#' ## website data test (printing with many decimals for 10 selected
#' ## questionnaires)
#' web <- c(1,2,4,5,11,27,28,31,37,39)
#' print(scores[web,], digits = 6)
#' ## SF12 Manual checks
#' print(unlist(lapply(scores, mean)), digits = 3)
#' print(unlist(lapply(scores, sd)), digits = 3)
#' print(lapply(scores, range), digits = 3)
#' ## Correlations
#' db <- cbind(sf12sample, scores)
#' var.order <- c(2:5,8,1,10,12,6,7,9,11)
#' cors <- cor(db)[var.order, 13:14]
#' print(cors, digits = 1)
#' ## Fine: reversed item have reverse sign correlation coefficients 
#' @export 
sf12 <- function( X = NULL ) {

  if((!(is.data.frame(X) | is.matrix(X))) | (ncol(X)!=12) )
    stop("X must be a data.frame (or matrix) with 12 columns")

  X <- as.data.frame(lapply(as.data.frame(X), as.integer))
  names(X) <- c("gh1", "pf2a", "pf2b", "rp3a", "rp3b", "re4a", "re4b", "bp5",
                "mh6a", "vt6b", "mh6c", "sf7" )
  
  ## *****************************************************************;
  ## ***               STEP 1: DATA CLEANING/REVERSE SCORING       ***;
  ## *****************************************************************;

  threept <- c("pf2a", "pf2b")
  fivept <- setdiff(names(X), threept)

  outRangeNA <- function(x, Min = 1L, Max) replace(x, x < Min | x > Max | is.null(x), NA)
  
  X[, threept] <- lapply(X[, threept], outRangeNA, Max = 3L)
  X[, fivept] <- lapply(X[, fivept], outRangeNA, Max = 5L)

  ghFunc <- function(i) {
    ghCalibrated <- list(1,2,3.4,4.4,5)
    return(ghCalibrated[[i]])
  }
  revFunc <- function(i) {
    bpSixabCalibrated <- list(5,4,3,2,1)
    return(bpSixabCalibrated[[i]])
  }

  X$ghc1 <- sapply(X$gh1, ghFunc)
  X$bpc5  <-  sapply(X$bp5, revFunc)
  X$mhc6a <- sapply(X$mh6a, revFunc)
  X$vtc6b <- sapply(X$vt6b, revFunc)

  ## *****************************************************************;
  ## *               STEP 2: CALCULATE RAW SCORES FROM               *
  ## *                       RECALIBRATED SCORES                     *
  ## *****************************************************************;
  
  pfRaw <- X$pf2a + X$pf2b
  rpRaw <- X$rp3a + X$rp3b
  bpRaw <- X$bpc5
  ghRaw <- X$ghc1
  vtRaw <- X$vtc6b
  sfRaw <- X$sf7
  reRaw <- X$re4a + X$re4b
  mhRaw <- X$mhc6a + X$mh6c

  ## *****************************************************************;
  ## *               STEP 3: SCALE RAW SCORES TO 0-100               * 
  ## *****************************************************************;

  scalePf <- function(rawScore) (rawScore - 2) / 6 * 100
  scaleRpReMh <- function(rawScore) (rawScore - 2) / 10 * 100
  scaleBpGhVtSf <- function(rawScore) (rawScore - 1) / 5 * 100


  pfScaled <- sapply(pfRaw, scalePf)
  rpScaled <- sapply(rpRaw, scaleRpReMh)
  bpScaled <- sapply(bpRaw, scaleBpGhVtSf)
  ghScaled <- sapply(ghRaw, scaleBpGhVtSf)
  vtScaled <- sapply(vtRaw, scaleBpGhVtSf)
  sfScaled <- sapply(sfRaw, scaleBpGhVtSf)
  reScaled <- sapply(reRaw, scaleRpReMh)
  mhScaled <- sapply(mhRaw, scaleRpReMh)


  ## *****************************************************************;
  ## *               STEP 4: STANDARDIZE SCALES WITH                 * 
  ## *                       Z-SCORE STANDARDIZATION                 * 
  ## *****************************************************************;

  pfStandardized <- (pfScaled - 81.18122) / 29.10558
  rpStandardized <- (rpScaled - 80.52856) / 27.13526
  bpStandardized <- (bpScaled - 81.74015) / 24.53019
  ghStandardized <- (ghScaled - 72.19795) / 23.19041
  vtStandardized <- (vtScaled - 55.59090) / 24.84380
  sfStandardized <- (sfScaled - 83.73973) / 24.75775
  reStandardized <- (reScaled - 86.41051) / 22.35543
  mhStandardized <- (mhScaled - 70.18217) / 20.50597

  ## *****************************************************************;
  ## *               STEP 5: NORM-BASED TRANSFORMATION OF            *
  ## *                       SCALE SCORES                            *
  ## *****************************************************************;

  nbsTransform <- function(s) {
      return(50 + 10 * s)
  }

  pfTransformed <- sapply(pfStandardized, nbsTransform)
  rpTransformed <- sapply(rpStandardized, nbsTransform)
  bpTransformed <- sapply(bpStandardized, nbsTransform)
  ghTransformed <- sapply(ghStandardized, nbsTransform)
  vtTransformed <- sapply(vtStandardized, nbsTransform)
  sfTransformed <- sapply(sfStandardized, nbsTransform)
  reTransformed <- sapply(reStandardized, nbsTransform)
  mhTransformed <- sapply(mhStandardized, nbsTransform)

  ## *****************************************************************;
  ## *               STEP 6: CALCULATE AGGREGATE SCORES              *
  ## *                       FOR PHYSICAL AND MENTAL                 *
  ## *****************************************************************;

  aggPhys <- 
    pfStandardized * 0.42402 +
    rpStandardized * 0.35119 +
    bpStandardized * 0.31754 +
    ghStandardized * 0.24954 +
    vtStandardized * 0.02877 +
    sfStandardized *-0.00753 +
    reStandardized *-0.19206 +
    mhStandardized *-0.22069

  aggMent <-
    pfStandardized *-0.22999 +
    rpStandardized *-0.12329 +
    bpStandardized *-0.09731 +
    ghStandardized *-0.01571 +
    vtStandardized * 0.23534 +
    sfStandardized * 0.26876 +
    reStandardized * 0.43407 +
    mhStandardized * 0.48581

  PCS_ <- sapply(aggPhys, nbsTransform)
  MCS_ <- sapply(aggMent, nbsTransform)

  result <- data.frame(pfNBS=pfTransformed,
                       rpNBS=rpTransformed,
                       bpNBS=bpTransformed,
                       ghNBS=ghTransformed,
                       vtNBS=vtTransformed,
                       sfNBS=sfTransformed,
                       reNBS=reTransformed,
                       mhNBS=mhTransformed,
                       PCS=PCS_,
                       MCS=MCS_)

  return(result)
}