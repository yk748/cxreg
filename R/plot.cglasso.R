#' plot heatmap from a "cglasso" object
#'
#' Produces plot of the estimated inverse spectral matrix for a fitted
#' \code{"cglasso"} object.
#'
#' A inverse spectral matrix profile plot is produced.
#'
#' @aliases plot.cglasso
#' @param X fitted \code{"cglasso"} model
#' @param index For which inverse spectral matrix profile is the plot to be drawn? Default is 1. 
#' The index must be provided within the length of the sequence of lambdas.
#' @param type Whether the plot is for real or imaginary part, or both, or in modulus (mod; absolute scale). Default is \code{mod}.
#' @param label If \code{TRUE}, label the axes with variable names.
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{yk748@cornell.edu}
#' @seealso \code{cglasso}
#'
#' @method plot cglasso
#' @export
plot.cglasso <- function(X, index, type=c("real","imaginary","mod","both"),label=FALSE) {
  

  ####################################################################
  # check for index for X
  if(!is.integer(index)){
    stop("The index must be integer.")
  }else if(1 > index | index > length(X)){
    stop("The index must be within the length of lambdas.")
  }
  
  if (is.null(type)){
    type <- "mod"
  }
  
  S <- X[[index]]
  p <- nrow(S)
  x_coords <- seq(0, 1, length.out = p + 1)  # breaks between columns
  y_coords <- seq(0, 1, length.out = p + 1)  # breaks between rows
  
  # variable labels
  if (label == TRUE){
    vnames <- colnames(S)
  }
  
  if (type != "both"){
    if (type == "real"){
      palette <- colorRampPalette(c("white", "blue"))
      S <- Re(S)
    }else if (type == "imaginary"){
      palette <- colorRampPalette(c("white", "red"))
      S <- Im(S)
    }else{
      palette <- colorRampPalette(c("white", "black"))
      S <- Mod(S)
    }
    
    # check type of X
    if (scale == "covariance"){
      z_scale <- c(min(S),max(S))
    }else{
      z_scale <- c(0,1)
    }
    
    S <- S[, nrow(S):1]
    image.plot(x = x_coords, y = y_coords,
               z = S,
               col = palette(100),
               zlim = z_scale,
               xaxt = "n",
               yaxt = "n",
               xlab = "",
               ylab = "",
               legend.shrink = 0.6,
               cex.main = 2)
    axis(1, at=1:p, labels=vars)
    axis(x, at=1:p, labels=vars)
    
  }else{
    palette_re <- colorRampPalette(c("white", "blue"))
    palette_im <- colorRampPalette(c("white", "red"))
    S_re <- Re(S); S_im <- Im(S)
    
    # check type of X
    if (scale == "covariance"){
      z_re <- c(min(S_re),max(S_re))
      z_im <- c(min(S_im),max(S_im))
    }else{
      z_re <- z_im <- c(0,1)
    }
    
    S_re <- S_re[, nrow(S_re):1]
    S_im <- S_im[, nrow(S_im):1]
    
    old_par <- par(mfrow = c(2, 1))
    image.plot(x = x_coords, y = y_coords,
               z = S_re,
               col = palette_re(100),
               zlim = z_re,
               xaxt = "n",
               yaxt = "n",
               xlab = "",
               ylab = "",
               legend.shrink = 0.6,
               cex.main = 2)
    axis(1, at=1:p, labels=vars)
    axis(x, at=1:p, labels=vars)
    image.plot(x = x_coords, y = y_coords,
               z = S_im,
               col = palette_im(100),
               zlim = z_re,
               xaxt = "n",
               yaxt = "n",
               xlab = "",
               ylab = "",
               legend.shrink = 0.6,
               cex.main = 2)
    axis(1, at=1:p, labels=vars)
    axis(x, at=1:p, labels=vars)
    par(old_par)
  }

}