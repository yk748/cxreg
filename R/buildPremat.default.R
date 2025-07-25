#' Build Prediction Matrix
#'
#' @export
buildPredmat.default <- function(outlist, lambda, x, foldid, alignment){


  predmat <- matrix(NA, nrow(x), length(lambda))
  nfolds <- max(foldid)
  nlams <- double(nfolds)
  nlambda <- length(lambda)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj <- outlist[[i]]
    # preds <- switch(alignment,
    #                 lambda=predict(fitobj, x[which, , drop = FALSE], s=lambda,,...),
    #                 fraction=predict(fitobj, x[which, , drop = FALSE],...) )
    preds <- x[which, , drop = FALSE] %*% fitobj$beta
    nlami <- min(ncol(preds),nlambda)
    predmat[which, seq(nlami)] <- preds[,seq(nlami)]
    if(nlami<nlambda){
      predmat[which,seq(from=nlami,to=nlambda)] <- preds[,nlami]
    }
  }
  rn <- rownames(x)
  sn <- paste("s",seq(0,length=nlambda),sep="")
  dimnames(predmat) <- list(rn,sn)

  predmat
}

#' Build Prediction Matrix
#'
#' @export
buildPredmat <- function(outlist, lambda, x, foldid, alignment){
  UseMethod("buildPredmat")
}

