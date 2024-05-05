# -------------------------------------------- #
# Soft threshold operator (typical lasso operator for a single entry)
soft_thrs_single <- function(z,thrs){
  
  ell_1 <- abs(z)
  s_thrs <- ( z*(1 - thrs/ell_1) )*ifelse(ell_1 >= thrs,1,0)
  
  return(s_thrs)
}

# -------------------------------------------- #
# Generalized soft threshold operator
gen_soft_thrs <- function(z,thrs){
  
  ell_2 <- norm(z,"2")
  s_thrs <- ( z*(1 - thrs/ell_2) )*ifelse(ell_2 >= thrs,1,0)
  
  return(s_thrs)
}

# -------------------------------------------- #