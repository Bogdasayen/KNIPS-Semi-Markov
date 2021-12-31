# Utility functions for the KNIPS hesim analysis

# Utility function to change NJR lower diagonal matrices into symmetric matrices
fill_in_symmetric <- function(s) {
  # Steps to convert it to a numeric matrix
  s <- matrix(as.numeric(as.matrix(s)), ncol = dim(s)[2], nrow = dim(s)[1])
  
  # Steps to fill in the upper triangular component with the lower triangle
  s.diag = diag(s)
  s[upper.tri(s,diag=T)] = 0
  s = s + t(s) + diag(s.diag)
  return(s)
}

# Utility function to permute order of variables in a covariance matrix
move_last_element_first <- function(s) {
  nrows <- dim(s)[1]
  ncols <- dim(s)[2]
  
  s.temp <- s
  s.temp[c(2:nrows), c(2:ncols)] <- s.temp[c(1:(nrows-1)), c(1:(ncols-1))]
  s.temp[1, ] <- s[nrows, c(ncols, c(1:(ncols - 1)))]
  s.temp[, 1] <- s[c(nrows, c(1:(nrows - 1))), ncols]
  return(s.temp)
}


read_rcs_covariance <- function(filename = NULL,
                                sheetname = NULL,
                              par_names = NULL) {
  rcs_covariance  <- as.matrix(read_excel(filename, sheet = sheetname))
  # Only look at the coefficients
  # Constant plus knot coefficients
  npars <- (dim(rcs_covariance)[2] - 1)/2 + 1
  rcs_covariance <- rcs_covariance[c(2:(1 + npars)), c(2:(1 + npars))]
  rcs_covariance <- fill_in_symmetric(rcs_covariance)
  # NJR outputs constant last but flexsurv/hesim needs it first
  # Need to move last element (constant) to front
  rcs_covariance <- move_last_element_first(rcs_covariance)
  
  # Create names if not supplied
  if(is.null(par_names)) par_names <- c("cons", paste0("rcs", 1:(npars - 1)))
  colnames(rcs_covariance) <- rownames(rcs_covariance) <-  par_names
  
  return(rcs_covariance)
}

