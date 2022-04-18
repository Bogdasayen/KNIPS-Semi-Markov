# Utility functions for the KNIPS hesim analysis

format_continuous_results <- function(x, n_digits = 2) {
  return(paste0(format(mean(x, na.rm = TRUE), digits = n_digits, nsmall = n_digits), ", (", 
                format(quantile(x, probs = 0.025, na.rm = TRUE), digits = n_digits, nsmall = n_digits), ", ",
                format(quantile(x, probs = 0.975, na.rm = TRUE), digits = n_digits, nsmall = n_digits), ")"))
}

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

summarise_results <- function(total_costs, total_qalys, intervention_names, reference_intervention) {
  results_table <- matrix(NA, nrow = length(intervention_names), ncol = 9) 
  colnames(results_table) <- c("Total costs", "Total QALYs", 
                               "Net benefit (£20k/QALY)", "Net benefit (£30k/QALY)",
                               "Incremental costs", "Incremental QALYs", 
                               "Incremental net benefit (£20k/QALY)", "Incremental net benefit (£30k/QALY)",
                               "ICER")
  rownames(results_table) <- intervention_names
  net_benefit_20k <- total_qalys * 20000 - total_costs
  net_benefit_30k <- total_qalys * 30000 - total_costs
  incremental_costs <- total_costs - total_costs[, which(intervention_names == reference_intervention)]
  incremental_qalys <- total_qalys - total_qalys[, which(intervention_names == reference_intervention)]
  # Summarise each intervention
  for(i in 1:length(intervention_names)) {
    results_table[i, "Total costs"] <- format_continuous_results(total_costs[, i])
    results_table[i, "Total QALYs"] <- format_continuous_results(total_qalys[, i])
    results_table[i, "Net benefit (£20k/QALY)"] <- format_continuous_results(net_benefit_20k[, i])
    results_table[i, "Net benefit (£30k/QALY)"] <- format_continuous_results(net_benefit_30k[, i])
    # Incremental results
    results_table[i, "Incremental costs"] <- format_continuous_results(incremental_costs[, i])
    results_table[i, "Incremental QALYs"] <- format_continuous_results(incremental_qalys[, i])
    results_table[i, "Incremental net benefit (£20k/QALY)"] <- format_continuous_results(net_benefit_20k[, i]  - net_benefit_20k[, which(intervention_names == reference_intervention)])
    results_table[i, "Incremental net benefit (£30k/QALY)"] <- format_continuous_results(net_benefit_30k[, i]  - net_benefit_30k[, which(intervention_names == reference_intervention)])
    # Calculate the ICER
    results_table[i, "ICER"] <- format(mean(incremental_costs[, i]) / mean(incremental_qalys[, i]), digits = 2, nsmall = 2)
    if(mean(incremental_costs[, i]) < 0 & mean(incremental_qalys[, i]) > 0) results_table[i, "ICER"] <- "Dominant"
    if(mean(incremental_costs[, i]) > 0 & mean(incremental_qalys[, i]) < 0) results_table[i, "ICER"] <- "Dominated"
  }
  # Make sure the reference intervention is the first row
  reference_row <- results_table[which(intervention_names == reference_intervention), ]
  results_table <- rbind(reference_row, results_table[ -which(intervention_names == reference_intervention), ])
  rownames(results_table)[1] <- reference_intervention
  results_table[1, c(5:9)] <- "-"
  
  return(results_table)
}

