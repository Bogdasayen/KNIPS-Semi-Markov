

generate_model_outputs <- function(hesim_dat, model_inputs) {
  
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
  
  
  ############################################################
  ## Simulate disease progression ############################
  ############################################################
  
  # Nested loop over implants and samples
  # Over implants allows independent spline models for each implant
  # Over samples is to allow 2nd revision to depend on time to 1st revision
  # disprog is a list of progressions for one implant with each element corresponding to one sample
  # disprog_combined_implant is a list of combined sampled progressions with each element correspdoning to one implant
  # disprog_combined is the combined progressions for all implants and samples
  # Only the last is kept due to memory limitations
  # Same for stateprobs, stateprobs_combined_implant and stateprobs_combined
  
  # List of disease progression and state occupancies, one for each implant
  disprog_combined_implant <- list()
  stateprobs_combined_implant <- list()
  
  # To test try implant_name <- "Cem CR_Fix Mono"
  system.time({
    for(implant_name in implant_names) {
      print(paste0("Implant ", which(implant_names == implant_name), "/", length(implant_names)))
      # Create a list of disease progression results for this implant
      disprog <- list()
      # And of state occupancies
      stateprobs <- list()
      
      # Simulate one at a time
      n_samples_temp <- 1
      # Simple function to convert vector to matrix with columns named constant ('cons')
      matrixv <- function(v, n = NULL){
        if (length(v) == 1) v <- rep(v, n_samples_temp) 
        m <- matrix(v)
        colnames(m) <- "cons"
        return(m)
      }
      
      for(i_sample in 1:n_samples) {
        #print(paste0("Sample ", i_sample))
        
        # Define the transition rate models
        
        # Function to generate random samples of underlying model parameters
        # Create a list of randomly sampled survival models for transitions
        transition_model_params <- params_surv_list(
          # 1. Post TKR to Post 1st revision
          params_surv(coefs = list(
            gamma0 = data.frame(int = rcs_first_revision[[implant_name]][i_sample, "cons"]),
            gamma1 = data.frame(int = rcs_first_revision[[implant_name]][i_sample, "rcs1"]),
            gamma2 = data.frame(int = rcs_first_revision[[implant_name]][i_sample, "rcs2"]),
            gamma3 = data.frame(int = rcs_first_revision[[implant_name]][i_sample, "rcs3"])),  
            aux = list(
              knots = ln_bhknots_first_revision[[implant_name]],
              scale = "log_hazard",
              timescale = "log"),
            dist = "survspline"),
          
          
          # 2. Post TKR to death
          params_surv(coefs = lapply(as.list(log(mr)), matrixv),
                      aux = list(time = c(0:time_horizon)),
                      dist = "pwexp"),
          
          # 3. Post 1st revision to post 2nd revision
          # Includes dependence on time from TKR to 1st revision
          params_surv(coefs = list(
            gamma0 = data.frame(int = rcs_second_revision[i_sample, "cons"],
                                yrs_to_1st_rev = rcs_second_revision[i_sample, "yrs_to_1st_rev"]),
            gamma1 = data.frame(int = rcs_second_revision[i_sample, "rcs1"]),
            gamma2 = data.frame(int = rcs_second_revision[i_sample, "rcs2"]),
            gamma3 = data.frame(int = rcs_second_revision[i_sample, "rcs3"])),  
            aux = list(
              knots = ln_bhknots_second_revision,
              scale = "log_hazard",
              timescale = "log"),
            dist = "survspline"),
          
          # 4. Post 1st revision to death
          params_surv(coefs = lapply(as.list(log(mr)), matrixv),
                      aux = list(time = c(0:time_horizon)),
                      dist = "pwexp"),
          
          # 5. Post 2nd revision to death
          params_surv(coefs = lapply(as.list(log(mr)), matrixv),
                      aux = list(time = c(0:time_horizon)),
                      dist = "pwexp"),
          
          # 6. Post TKR to surgery death
          params_surv(coefs = list(rate = log(tkr_mortality_90d[i_sample, 1]),
                                   as.matrix(rep(-100, n_samples_temp))),
                      aux = list(time = tkr_mortality_times_90d),
                      dist = "pwexp"),
          
          # 7. Post 1st revision to surgery death
          params_surv(coefs = list(rate = log(revision_mortality_90d[i_sample, 1]),
                                   as.matrix(rep(-100, n_samples_temp))),
                      aux = list(time = revision_mortality_times_90d),
                      dist = "pwexp"),
          
          # 8. Post 2nd revision to surgery death
          # Depends on initial 90d period and ongoing risk of higher revision
          params_surv(coefs = list(rerevision_mortality_90d[i_sample, 1],
                                   rerevision_mortality_90d[i_sample, 2]),
                      aux = list(time = rerevision_mortality_times_90d),
                      dist = "pwexp")
          
        )
        
        # Now simulate the outcomes for all patients for this sample and this implant
        # Need temporary hesim data object with only one strategy
        hesim_dat_temp <- hesim_data(strategies = strategies[strategies$strategy_name == implant_name],
                                     patients = patients, 
                                     states = states)
        
        transition_model_data <- expand(hesim_dat_temp, by = c("strategies", "patients"))
        
        # Add numerical covariates
        # Need to add constants to match the spline and piecewise constant death parameters
        transition_model_data[, cons := 1]
        transition_model_data[, int := 1]
        transition_model_data[, x1 := 1]
        transition_model_data[, yrs_to_1st_rev := 0]
        
        # Transition model
        transition_model <- create_IndivCtstmTrans(transition_model_params, 
                                           input_data = transition_model_data,
                                           trans_mat = tmat,
                                           clock = "mix",
                                           start_age = patients$age,
                                           reset_states = c(1, 2))
        
        
        set.seed(2243534)
        disprog[[i_sample]] <- transition_model$sim_disease(max_t = time_horizon, max_age = (starting_age + time_horizon))
        
        # Just use to test if state probabilities are making sense
        #x <-  transition_model$sim_stateprobs(t = 0:(time_horizon))
        # with(x, prob[strategy_id ==1 & sample == 1 & t==39])
        
        # Update the data using results of this initial simulation
        # Use the time to transition from state 1 to state 2
        # Not all patients on all transitions actually move from 1 to 2
        disprog_temp <- disprog[[i_sample]][from == 1 & to == 2, ]
        
        # Using slow loop but might be better way
        for(i_transition in 1:dim(disprog_temp)[1]) {
          transition_model_data[transition_model_data$patient_id == as.numeric(disprog_temp[i_transition, "patient_id"]) &
                          transition_model_data$strategy_id == as.numeric(disprog_temp[i_transition, "strategy_id"]), "yrs_to_1st_rev"] <-
            disprog_temp[i_transition, "time_stop"]
          
        }
        
        # Re-run simulation with the (correct) time to first revision
        # Transition model
        transition_model <- create_IndivCtstmTrans(transition_model_params, 
                                           input_data = transition_model_data,
                                           trans_mat = tmat,
                                           clock = "mix",
                                           start_age = patients$age,
                                           start_state = 1,
                                           reset_states = c(1, 2))
        set.seed(2243534)
        disprog[[i_sample]] <- transition_model$sim_disease(max_t = time_horizon, max_age = (starting_age + time_horizon))
        # Set the sample number
        disprog[[i_sample]]$sample <- i_sample
        disprog[[i_sample]]$strategy_id <- which(implant_names == implant_name)
        
        # Also simulate average state probabilities
        stateprobs[[i_sample]] <- transition_model$sim_stateprobs(t = 0:(time_horizon))
        stateprobs[[i_sample]]$sample <- i_sample
        stateprobs[[i_sample]]$strategy_id <- which(implant_names == implant_name)
        
        # Check that state probs make sense
        # with(stateprobs[[i_sample]], prob[strategy_id ==1 & sample == 1 & t==20])
      } # End loop over samples
      
      
      
      # Create disease progression combined for this implant
      disprog_combined_implant[[implant_name]] <-rbindlist(disprog)
      stateprobs_combined_implant[[implant_name]] <- rbindlist(stateprobs)
      
      # Correctly set the size attributes
      # Check if using attributes(disprog_combined)
      attributes_size <- c(n_samples, n_strategies, n_patients, n_states + 2)
      names(attributes_size) <- c("n_samples", "n_strategies", "n_patients", "n_states")
      setattr(disprog_combined_implant[[implant_name]], "size", attributes_size)
      # Correctly set the absorbing attribute
      attributes_absorbing <- c(4, 5)
      names(attributes_absorbing) <- c("Surgery Death", "Death")
      setattr(disprog_combined_implant[[implant_name]], "absorbing", attributes_absorbing)
      
    } # End loop over implants
  }) # End system time
  
  # Create disprog_combined
  disprog_combined <-rbindlist(disprog_combined_implant)
  stateprobs_combined <- rbindlist(stateprobs_combined_implant)
  
  # Check probabilities add to one
  # with(stateprobs[[1]], prob[sample == 1 & strategy_id == 12 & t == 15])
  # with(stateprobs_combined_implant[[12]], prob[sample == 1 & strategy_id == 12 & t == 15])
  # with(stateprobs_combined, prob[sample == 1 & strategy_id == 12 & t == 15])
  
  
  # Correctly set the size attributes
  # Check if using attributes(disprog_combined)
  attributes_size <- c(n_samples, n_strategies, n_patients, n_states + 2)
  names(attributes_size) <- c("n_samples", "n_strategies", "n_patients", "n_states")
  setattr(disprog_combined, "size", attributes_size)
  # Correctly set the absorbing attribute
  attributes_absorbing <- c(4, 5)
  names(attributes_absorbing) <- c("Surgery Death", "Death")
  setattr(disprog_combined, "absorbing", attributes_absorbing)
  
  ############################################################
  ## Combined economic model #################################
  ############################################################
  
  
  cost_utility_models <- generate_cost_utility_models(n_samples = n_samples,
                                                      hesim_dat = hesim_dat,
                                                      model_inputs = model_inputs,
                                                      disprog_combined = disprog_combined)
  
  # Economic model combining everything
  economic_model <- IndivCtstm$new(trans_model = transition_model,
                            utility_model = cost_utility_models$utility_model,
                            cost_models = cost_utility_models$cost_models)
  
  # Use the combined disease progression results
  economic_model$disprog_ <- disprog_combined
  
  # Look at simulations to ensure the combination was correct
  head(economic_model$disprog_)
  
  # And the combined state occupancy probabilities
  economic_model$stateprobs_ <- stateprobs_combined
  
  economic_model$stateprobs_
  # Average probability of occupying Post TKR on strategy 1
  with(economic_model$stateprobs_, mean(prob[strategy_id == 1 & state_id == 1]))
  # Average in post 1st revision
  with(economic_model$stateprobs_, mean(prob[strategy_id == 1 & state_id == 2]))
  # Average in post 2nd revision
  with(economic_model$stateprobs_, mean(prob[strategy_id == 1 & state_id == 3]))
  # Average probability of occupying surgery death
  with(economic_model$stateprobs_, mean(prob[strategy_id == 1 & state_id == 4]))
  # Average probability of occupying death
  with(economic_model$stateprobs_, mean(prob[strategy_id == 1 & state_id == 5]))
  
  # Simulate QALYs and costs with discount rates
  # 3.5% discount rate
  economic_model$sim_qalys(dr = discount_rate)
  economic_model$sim_costs(dr = discount_rate)
  
  model_outputs <- list("economic_model" = economic_model)
  
}