#' soft_thrs_single
#'
#' Soft threshold operator (typical lasso operator for a single entry)
#' @param z input1
#' @param thrs input2
#' @return output
#' @examples
#' example <- soft_thrs_single(z,thrs)
#'
#' @export
soft_thrs_single <- function(z,thrs){

  ell_1 <- abs(z)
  s_thrs <- ( z*(1 - thrs/ell_1) )*ifelse(ell_1 >= thrs,1,0)

  return(s_thrs)
}



#' soft_thrs_single
#'
#' Generalized soft threshold operator
#' @param z input1
#' @param thrs input2
#' @return output
#' @examples
#' example <- gen_soft_thrs(z,thrs)
#'
#' @export
gen_soft_thrs <- function(z,thrs){

  ell_2 <- norm(z,"2")
  s_thrs <- ( z*(1 - thrs/ell_2) )*ifelse(ell_2 >= thrs,1,0)

  return(s_thrs)
}

# -------------------------------------------- #
