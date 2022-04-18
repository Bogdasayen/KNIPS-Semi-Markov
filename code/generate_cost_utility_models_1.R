

generate_cost_utility_models <- function(n_samples, 
                                         hesim_dat,
                                         model_inputs,
                                         disprog_combined) {
  
  # Extract requisite random and fixed parameters from model_inputs
  # This is to ensure input_parameters aligns with format required for BCEA-EVPPI
  
  # Fixed parameters
  ln_bhknots_first_revision <- model_inputs$ln_bhknots_first_revision
  ln_bhknots_second_revision <- model_inputs$ln_bhknots_second_revision
  mr <- model_inputs$mr
  tkr_mortality_times_90d <- model_inputs$tkr_mortality_times_90d
  revision_mortality_times_90d <- model_inputs$revision_mortality_times_90d
  rerevision_mortality_times_90d <- model_inputs$rerevision_mortality_times_90d
  
  tkr_surgery_cost <- model_inputs$tkr_surgery_cost
  revision_cost <- model_inputs$revision_cost
  
  # Random parameters
  input_parameters <- model_inputs$inputs_parameters
  rcs_second_revision <- input_parameters[, grep("rcs_second_revision", colnames(input_parameters))]
  colnames(rcs_second_revision) <- gsub("rcs_second_revision_", "", colnames(rcs_second_revision))
  
  probability_higher_revision <- input_parameters[, "probability_higher_revision"]
  
  rcs_first_revision <- list()
  for(implant_name in implant_names) {
    temp <- input_parameters[, grep(paste0("rcs_first_revision_", implant_name), colnames(input_parameters))]
    
    colnames(temp) <- gsub(paste0("rcs_first_revision_", implant_name, "_"), "", colnames(temp))
    rcs_first_revision[[implant_name]] <- temp
  }
  
  tkr_mortality_90d <- input_parameters[, grep("tkr_mortality_90d_", colnames(input_parameters))]
  revision_mortality_90d <- input_parameters[, grep("revision_mortality_90d_", colnames(input_parameters))]
  rerevision_mortality_90d <- input_parameters[, grep("rerevision_mortality_90d_", colnames(input_parameters))]
  
  revision_disutility <- input_parameters[, "revision_disutility"]
  utility_post_tkr <- input_parameters[, "utility_post_tkr"]
  utility_post_1st_rev <- input_parameters[, "utility_post_1st_rev"]
  utility_post_2nd_rev <- input_parameters[, "utility_post_2nd_rev"]
  
  implant_costs <- input_parameters[, grep("implant_cost_", colnames(input_parameters))]
  colnames(implant_costs) <- gsub("implant_cost_", "", colnames(implant_costs))
  implant_costs <- t(implant_costs)
  
  
  # Need probabilities of 1st and 2nd revision in each 1-year time interval
  # Probability 1st revision depends on strategy, sample, and time (not patient)
  probability_1st_revision <- data.table(strategy_id = rep(strategies$strategy_id, each = (n_samples * n_time_intervals)),
                                         sample = rep(rep(1:n_samples, times = n_strategies), each = n_time_intervals),
                                         time_start = rep(0:(n_time_intervals - 1), times = (n_strategies * n_samples)),
                                         value = rep(0, times = n_strategies * n_samples * n_time_intervals))
  
  
  # Probability 2nd revision depends on strategy, patient, sample, and time 
  probability_2nd_revision <- data.table(strategy_id = rep(strategies$strategy_id, each = (n_patients * n_samples * n_time_intervals)),
                                         patient_id = rep(rep(c(1:n_patients), times = n_strategies), each = (n_samples * n_time_intervals)),
                                         sample = rep(rep(1:n_samples, times = n_strategies * n_patients), each = n_time_intervals),
                                         time_start = rep(0:(n_time_intervals - 1), times = (n_strategies * n_patients * n_samples)),
                                         value = rep(0, times = n_strategies * n_patients * n_samples * n_time_intervals))
  
  # For efficiency, create a single matrix with spline parameters
  rcs_first_revision_temp <- rbindlist(rcs_first_revision)
  # And a matrix with one set of knots for each sample
  ln_bhknots_first_revision_temp <- rbindlist(lapply(ln_bhknots_first_revision, as.data.frame))
  ln_bhknots_first_revision_temp <- ln_bhknots_first_revision_temp[rep(1:12, each = n_samples), ]
  
  # Intervals are 1 year if time_horizon is equal to n_time_intervals
  for(time_start_ in 0:(n_time_intervals - 1)) {
    # Varies by strategy and sample
    probability_1st_revision[probability_1st_revision$time_start == time_start_, "value"] <- 
      exp(-Hsurvspline(
        x = time_start_,
        gamma = as.matrix(rcs_first_revision_temp),
        knots = as.matrix(ln_bhknots_first_revision_temp))) - 
      exp(-Hsurvspline(
        x = time_start_ + time_horizon/n_time_intervals,
        gamma = as.matrix(rcs_first_revision_temp),
        knots = as.matrix(ln_bhknots_first_revision_temp)
      ))
  }
  
  # Depends on patient, strategy and sample
  # Using 100 as a placeholder for rows that have no 1st revision
  # These placeholders aren't used in the analysis
  yrs_to_1st_rev_temp <- data.table(strategy_id = rep(strategies$strategy_id, each = (n_patients * n_samples)),
                                    patient_id = rep(rep(c(1:n_patients), times = n_strategies), each = n_samples),
                                    sample = rep(1:n_samples, times = n_strategies * n_patients),
                                    time_stop = rep(100, times = n_strategies * n_patients * n_samples))
  
  # Rows with a 1st revision
  disprog_temp <- disprog_combined[disprog_combined$from == 1 & disprog_combined$to == 2, c("strategy_id", "patient_id", "sample", "time_stop")] 
  
  # Fill in the gaps (i.e. patients/strategies for which transition doesn't occur, and don't affect model results)
  # Set the rows of from placeholders to actual 1st revision times
  setkey(yrs_to_1st_rev_temp, strategy_id, patient_id, sample)
  for(i in 1:dim(disprog_temp)[1]) {
    yrs_to_1st_rev_temp[.(disprog_temp[i, "strategy_id"], 
                          disprog_temp[i, "patient_id"], 
                          disprog_temp[i, "sample"]), "time_stop"] <- 
      disprog_temp[i, "time_stop"]
  }
  
  
  # Probabilities of 2nd revision in Post 1st revision state
  # Intervals are 1 year if time_horizon equals n_time_horizon
  
  # Spline parameters and the covariate effect are separated for Hsurvspline
  gamma_temp <- rcs_second_revision[, -2]
  beta_temp <- rcs_second_revision[, 2]
  
  
  
  # Create duplicates for each patient and strategy
  # These only vary by sample
  gamma_temp <- gamma_temp[rep(c(1:n_samples), n_patients * n_strategies), ]
  beta_temp <- rep(beta_temp, n_patients * n_strategies)
  
  # Ensure it is a matrix
  gamma_temp <- as.matrix(gamma_temp)
  
  for(time_start_ in 0:(n_time_intervals - 1)) {
    print(paste0("Time interval ", time_start_, "/", n_time_intervals))
    
    # Hsurvspline expects x and X to vectors of times and covariate values, respectively
    # gamma is a matrix with one row for each sample
    # beta is a vector with on element for each sample
    # In code below the x and X are fixed and effect*covariate is calculated before passing to Hsurvspline
    # Tested element by element to ensure gives correct answer
    probability_2nd_revision[probability_2nd_revision$time_start == time_start_, "value"] <- 
      exp(-Hsurvspline(
        x = time_start_,
        knots = ln_bhknots_second_revision,
        gamma = gamma_temp,
        beta = beta_temp * yrs_to_1st_rev_temp$time_stop,
        X = 1)) - 
      exp(-Hsurvspline(
        x = time_start_ + time_horizon/n_time_intervals,
        knots = ln_bhknots_second_revision,
        gamma = gamma_temp,
        beta = beta_temp * yrs_to_1st_rev_temp$time_stop, 
        X = 1)) 
    
  } # End loop over times
  
  
  
  
  # Utilities depend on the probabilities of revision
  # Default value is zero so utilities in Surgery Death are zero
  # This may not work if n_patients or n_samples becomes too big
  utility_tbl <- stateval_tbl(
    data.table(strategy_id = rep(strategies$strategy_id, each = (n_patients * n_states * n_samples * n_time_intervals)),
               patient_id = rep(rep(1:n_patients, times = n_strategies), each = (n_states * n_samples * n_time_intervals)),
               state_id = rep(rep(states$state_id, times = n_strategies * n_patients), each = (n_samples * n_time_intervals)),
               sample = rep(rep(1:n_samples, times = n_strategies * n_patients * n_states), each = n_time_intervals),
               time_start = rep(0:(n_time_intervals - 1), times = (n_strategies * n_patients * n_states * n_samples)),
               value = rep(0, times = n_strategies * n_patients * n_states * n_samples * n_time_intervals)),
    dist = "custom"
  )
  
  
  # Utilities in Post TKR are dependent on time (clock-forward as all start in Post TKR)
  for(time_start_ in 0:(n_time_intervals - 1)) {
    utility_tbl[utility_tbl$state_id == 1 & utility_tbl$time_start == time_start_, "value"] <- rep(utility_post_tkr, n_strategies * n_patients) +# general utility
      # Disutility times probability of revision during this time interval
      rep(revision_disutility, n_strategies * n_patients) *
      rep(unlist(probability_1st_revision[probability_1st_revision$time_start == time_start_, "value"]), times = n_patients)
  }
  
  
  
  
  # Utilities in Post 1st revision are dependent on time 
  # Note that this is clock reset
  for(time_start_ in 0:(n_time_intervals - 1)) {
    utility_tbl[utility_tbl$state_id == 2 & utility_tbl$time_start == time_start_, "value"] <- rep(utility_post_1st_rev, n_strategies * n_patients) +# general utility
      # Disutility times probability of revision during this time intervals
      rep(revision_disutility, n_strategies * n_patients) * 
      probability_2nd_revision[probability_2nd_revision$time_start == time_start_, "value"]
  }
  
  # Utility post 2nd revision 
  # Don't depend on time, time to second revision, or strategy
  # The subtraction may not be needed if utility_post_2nd_rev already includes consequence of higher revisions
  utility_tbl[utility_tbl$state_id == 3, "value"] <- 
    rep(utility_post_2nd_rev + probability_higher_revision * revision_disutility, 
        each = n_strategies * n_patients * n_time_intervals) 
  
  
  
  # Check that all strategies currently have the same utility in Post TKR
  #utility_tbl[utility_tbl$state_id == 1 & utility_tbl$time_start == 20 &
  #             utility_tbl$sample == 1, ]
  # And in Post 1st revision
  #utility_tbl[utility_tbl$state_id == 2 & utility_tbl$time_start == 20 &
  #     utility_tbl$sample == 1, ]
  
  
  head(utility_tbl)
  
  
  # Implant costs only in "Post TKR" state and depends on strategy
  implantcost_tbl <- stateval_tbl(
    data.table(strategy_id = rep(strategies$strategy_id, each = n_states * n_samples),
               state_id = rep(rep(states$state_id, times = n_strategies), each = n_samples),
               sample = rep(1:n_samples, times = n_states * n_strategies),
               value = 0),
    dist = "custom"
  )
  
  for(i_implant in 1:n_strategies) {
    implantcost_tbl[implantcost_tbl$strategy_id == i_implant &
                      implantcost_tbl$state_id == 1, "value"] = tkr_surgery_cost + implant_costs[i_implant, ]
  }
  
  
  # Revision costs in all states and depend on strategy, sample and patient
  # May not work if too many patients or samples
  medcost_tbl <- stateval_tbl(
    data.table(strategy_id = rep(strategies$strategy_id, each = (n_patients * n_states * n_samples * n_time_intervals)),
               patient_id = rep(rep(1:n_patients, times = n_strategies), each = (n_states * n_samples * n_time_intervals)),
               state_id = rep(rep(states$state_id, times = n_strategies * n_patients), each = (n_samples * n_time_intervals)),
               sample = rep(rep(1:n_samples, times = n_strategies * n_patients * n_states), each = n_time_intervals),
               time_start = rep(0:(n_time_intervals - 1), times = (n_strategies * n_patients * n_states * n_samples)),
               value = rep(0, times = n_strategies * n_patients * n_states * n_samples * n_time_intervals)),
    dist = "custom"
  )
  
  # Costs in Post TKR are dependent on time (clock-forward as all start in Post TKR)
  for(time_start_ in 0:(n_time_intervals - 1)) {
    medcost_tbl[medcost_tbl$state_id == 1 & medcost_tbl$time_start == time_start_, "value"] <-  
      # Surgery cost times probability of revision during this time intervals
      + rep(revision_cost, n_strategies * n_patients) *
      rep(unlist(probability_1st_revision[probability_1st_revision$time_start == time_start_, "value"]), times = n_patients)
    
  }
  
  # Costs in Post 1st revision are dependent on time 
  # Note this is clock reset
  for(time_start_ in 0:(n_time_intervals - 1)) {
    medcost_tbl[medcost_tbl$state_id == 2 & medcost_tbl$time_start == time_start_, "value"] <- 
      # Cost times probability of revision during this time interval
      + rep(revision_cost, n_strategies * n_patients) * 
      probability_2nd_revision[probability_2nd_revision$time_start == time_start_, "value"]
  }
  
  # Cost post 2nd revision the same for all strategies, patients and times
  # It is revision cost times probability of higher revision
  medcost_tbl[medcost_tbl$state_id == 3, "value"] <- rep(revision_cost * probability_higher_revision, each = n_strategies * n_patients * n_time_intervals) 
  
  
  ############################################################
  ## Create cost and utility models
  ############################################################
  
  # Utility and cost models
  # Utility
  utility_model <- create_StateVals(utility_tbl, n = n_samples, 
                                 time_reset = TRUE) # Not totally sure this argument is needed
  
  # Costs
  # The 'starting' option means costs are only applied when patients enter the state
  implant_cost_model <- create_StateVals(implantcost_tbl, n = n_samples,
                                     method = "starting",
                                     hesim_data = hesim_dat) # Expand by patients
  
  medical_cost_model <- create_StateVals(medcost_tbl, n = n_samples,
                                 time_reset = TRUE) # Not totally sure this argument is needed
  
  cost_models <- list(Drug = implant_cost_model,
                   Medical = medical_cost_model)
  
  return(list("utility_model" = utility_model, 
              "cost_models" = cost_models))
  
}
