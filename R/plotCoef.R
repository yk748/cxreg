plotCoef <- function(beta,norm,lambda,df,dev,label=FALSE,
                     xvar=c("norm","lambda"),xlab=iname,...){

  ##beta should be in "dgCMatrix" format
  which <- nonzeroCoef(beta)
  nwhich <- length(which)
  switch(nwhich+1,#we add one to make switch work
         "0"={warning("No plot produced since all coefficients zero")
              return()
         },
         "1"=warning("1 or less nonzero coefficients; classo plot is not meaningful")
  )
  beta <- as.matrix(beta[which,,drop=FALSE])

  # ------------------------------------------------ #
  xvar<-match.arg(xvar)
  switch(xvar,
         "norm"={
           index=if(missing(norm))apply(abs(beta),2,sum)else norm
           iname="L1 Norm"
           approx.f=1
         },
         "lambda"={
           index=log(lambda)
           iname="Log Lambda"
           approx.f=0
         }
         # "dev"= {
         #   index=dev
         #   iname="Fraction Deviance Explained"
         #   approx.f=1
         # }
  )

  # ------------------------------------------------ #
  dotlist <- list(...)
  type <- dotlist$type
  if(is.null(type)){
    par(mfrow=c(2,1),mar=c(3,4,1,2))
    matplot(index,t(Re(beta)),lty=1,xlab=xlab,ylab="Coefficients (Re)",type="l",...)
    matplot(index,t(Im(beta)),lty=1,xlab=xlab,ylab="Coefficients (Im)",type="l",...)
  } else {
    par(mfrow=c(2,1),mar=c(3,4,1,2))
    matplot(index,t(Re(beta)),lty=1,xlab=xlab,ylab="Coefficients (Re)",...)
    matplot(index,t(Im(beta)),lty=1,xlab=xlab,ylab="Coefficients (Im)",...)
  }
  par()
  
  atdf <- pretty(index)

  # ------------------------------------------------ #
  ### compute df by interpolating to df at next smaller lambda
  prettydf <- approx(x=index,y=df,xout=atdf,rule=2,method="constant",f=approx.f)$y

  # prettydf=ceiling(approx(x=index,y=df,xout=atdf,rule=2)$y)
  axis(3,at=atdf,labels=prettydf,tcl=NA)
  if(label){
    nnz <- length(which)
    xpos <- max(index)
    pos <- 4
    if(xvar=="lambda"){
      xpos <- min(index)
      pos <- 2
    }
    xpos <- rep(xpos,nnz)
    ypos <- abs(beta[,ncol(beta)])
    text(xpos,ypos,paste(which),cex=.5,pos=pos)
  }

}
