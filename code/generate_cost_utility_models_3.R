 
generate_cost_utility_models <- function(n_samples, 
                                         hesim_dat,
                                         model_inputs) {
  
  # Extract requisite random and fixed parameters from model_inputs
  # This is to ensure input_parameters aligns with format required for BCEA-EVPPI
  
  tkr_surgery_cost <- model_inputs$tkr_surgery_cost
  revision_cost <- model_inputs$revision_cost
  
  # Random parameters
  input_parameters <- model_inputs$inputs_parameters
  probability_higher_revision <- input_parameters[, "probability_higher_revision"]
  
  revision_disutility <- input_parameters[, "revision_disutility"]
  utility_post_tkr <- input_parameters[, "utility_post_tkr"]
  utility_post_1st_rev <- input_parameters[, "utility_post_1st_rev"]
  utility_post_2nd_rev <- input_parameters[, "utility_post_2nd_rev"]
  
  implant_costs <- input_parameters[, grep("implant_cost_", colnames(input_parameters))]
  colnames(implant_costs) <- gsub("implant_cost_", "", colnames(implant_costs))
  implant_costs <- t(implant_costs)
  
  
  
  utility_tbl <- stateval_tbl(
    data.table(sample = rep(1:n_samples, each = n_states),
               state_id = rep(states$state_id, times = n_samples),
               value = rep(0, times = n_states * n_samples)),
    dist = "custom"
  )
  
  utility_tbl[state_id == 1, "value"] <- utility_post_tkr
  utility_tbl[state_id == 2, "value"] <- utility_post_1st_rev
  # Utility post 2nd revision 
  # Need to subtract disutilities of higher revision multiplied by (annual) probability of having revision
  # Don't depend on time, time to second revision, or strategy
  utility_tbl[state_id == 3, "value"] <- utility_post_2nd_rev + 
    probability_higher_revision * revision_disutility
            
  
  
  
  # Implant costs only in "Post TKR" state and depends on strategy
  implantcost_tbl <- stateval_tbl(
    data.table(strategy_id = rep(strategies$strategy_id, each = n_states * n_samples),
               sample = rep(rep(1:n_samples, each = n_states), times = n_strategies),
               state_id = rep(states$state_id, times = n_strategies * n_samples),
               value = 0),
    dist = "custom"
  )
  
  # Only non-zero for first state
  for(i_implant in 1:n_strategies) {
    implantcost_tbl[implantcost_tbl$strategy_id == i_implant &
                      implantcost_tbl$state_id == 1, "value"] = tkr_surgery_cost + implant_costs[i_implant, ]
  }
  
  
  # Medical costs
  medcost_tbl <- stateval_tbl(
    data.table(sample = rep(1:n_samples, each = n_states),
               state_id = rep(states$state_id, times = n_samples),
               value = rep(0, times = n_states * n_samples)),
    dist = "custom"
  )
  
  # No cost in post-TKR
  medcost_tbl[state_id == 1, "value"] <- 0
  # No cost in post 1st revision
  medcost_tbl[state_id == 2, "value"] <- 0
  # Cost post 2nd revision 
  # Cost of a higher revision multiplied by (annual) probability of having revision
  # Don't depend on time, time to second revision, or strategy
  medcost_tbl[state_id == 3, "value"] <-  probability_higher_revision * revision_cost
  

  ############################################################
  ## Create cost and utility models
  ############################################################
  
  # Utility and cost models
  # Utility
  utility_model <- create_StateVals(utility_tbl, n = n_samples, 
                                 time_reset = TRUE,
                                 hesim_data = hesim_dat) # Not totally sure this argument is needed
  
  # Costs
  # The 'starting' option means costs are only applied when patients enter the state
  implant_cost_model <- create_StateVals(implantcost_tbl, n = n_samples,
                                     method = "starting",
                                     hesim_data = hesim_dat) # Expand by patients
  
  medical_cost_model <- create_StateVals(medcost_tbl, n = n_samples,
                                 time_reset = TRUE,
                                 hesim_data = hesim_dat) 
  
  cost_models <- list(Drug = implant_cost_model,
                   Medical = medical_cost_model)
  
  return(list("utility_model" = utility_model, 
              "cost_models" = cost_models))
  
}
