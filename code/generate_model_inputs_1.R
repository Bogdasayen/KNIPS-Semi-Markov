# Generate an inputs list for the economic model
# This includes a matrix of random sampled parameters for use in BCEA-EVPPI calculations
# Fixed values are not included in the input_parameters matrix but as separate elements in model_inputs

generate_model_inputs <- function(n_samples, age_range, sample_gender, use_dummy_data = FALSE) {
  
  
  # Read in UK lifetables
  lifetables <- read_excel(paste0(data_directory, "/KNIPS Main input data.xlsx"), sheet = "uk_lifetables")
  
  
  ############################################################
  ## Transition rate parameters ##############################
  ############################################################
  
  
  # Transition 1
  # From Post TKR to Post 1st revision
  
  if(use_dummy_data) {
    first_revision_filename <- paste0("first_revision_test.xlsx")
  }  else {
    if(age_range[1] != 0 & !is.infinite(age_range[2])) {
      first_revision_filename <- paste0("/", age_range[1], "-", age_range[2], "_", sample_gender, "_first_revision.xlsx")  
    } 
    if(age_range[1] == 0) {
      # Less than age_range[2]
      first_revision_filename <- paste0("/", "l", age_range[2], "_", sample_gender, "_first_revision.xlsx")  
    }
    if(is.infinite(age_range[2])) {
      # Above age_range[1]
      first_revision_filename <- paste0("/", "a", age_range[1], "_", sample_gender, "_first_revision.xlsx")  
    }
  }
  
  # List of sampled spline parameters and (fixed) knot locations, one for each implant
  rcs_first_revision <- list()
  ln_bhknots_first_revision <- list()
  
  # For test try implant_name <- "Cem CR_Fix Mono"
  for(implant_name in implant_names) {
    rcs_first_revsision_mean_raw <- as.data.frame(read_excel(paste0(data_directory,first_revision_filename), sheet = "rcs_first_revision_mean"))
    ln_bhknots_first_revision_raw = as.data.frame(read_excel(paste0(data_directory,first_revision_filename), sheet = "ln_bhknots_first_revision"))
    rownames(rcs_first_revsision_mean_raw) <- rownames(ln_bhknots_first_revision_raw ) <- rcs_first_revsision_mean_raw[, 1]
    
    rcs_first_revision_mean <- as.matrix(rcs_first_revsision_mean_raw[implant_name, -1])
    ln_bhknots_first_revision[[implant_name]] <- as.matrix(ln_bhknots_first_revision_raw[implant_name, -1])
    
    # This is the covariance matrix provided by Linda from NJR in log_2_1
    rcs_first_revision_covariance <- read_rcs_covariance(filename = paste0(data_directory,first_revision_filename),
                                                         sheetname = paste0(implant_name, "_cov"),
                                                         par_names = names(rcs_first_revision_mean))
    
    
    rcs_first_revision[[implant_name]] <- as.data.frame(mvrnorm(n_samples, mu = rcs_first_revision_mean,
                                                                Sigma = rcs_first_revision_covariance))
    
    # Ensure it is a matrix
    if(n_samples == 1) {
      rcs_first_revision[[implant_name]] <- t(as.matrix(rcs_first_revision[[implant_name]], nrow = 1))
    }
  }
  
  
  # Rates of 3rd or higher revision
  # Used for surgery mortality after 2nd revision (transition 8) and in costs and utilities
  # Inverse variance meta-analysis of logrates of 3rd to 8th revision from NJR
  other_rates_raw <- as.matrix(read_excel(paste0(data_directory, "/KNIPS Main input data.xlsx"), sheet = "other_rates"))
  rownames(other_rates_raw) <- other_rates_raw[, "parameter"]
  
  lograte_higher_revision_mean <-  as.numeric(other_rates_raw["lograte_higher_revision_mean", "value"])
  lograte_higher_revision_se <-  as.numeric(other_rates_raw["lograte_higher_revision_se", "value"])
  
  # Assumed not to depend on time, time to 2nd revision, or implant
  # Need annualised probabilities for multiplying by 90d mortality and in costs and utilities in "Post 2nd revision"
  probability_higher_revision <- 1 - exp(- exp(rnorm(n = n_samples,
                                                     mean = lograte_higher_revision_mean,
                                                     sd = lograte_higher_revision_se)))
  
  
  
  # Transition 2
  # Post TKR to death 
  # Currenly background only but need to include KNIPS data
  # Annualised mortality probability
  mr <- as.matrix(lifetables[starting_age:(starting_age + time_horizon), paste0(sample_gender, "s")])
  mr_times <- c(0:time_horizon)
  
  # Transition 3
  # Post 1st revision to Post 2nd revision
  # These do not depend on age or gender  but do depend on time to 1st revision
  
  
  # NJR Spline model with 3 knots from "Re-revisions for main analysis_corrected.docx"
  # Parameters
  rcs_second_revision_mean <- as.matrix(read_excel(paste0(data_directory,"/KNIPS Main input data.xlsx"), sheet = "rcs_second_revision_mean"))
  # Knots
  ln_bhknots_second_revision = as.matrix(read_excel(paste0(data_directory,"/KNIPS Main input data.xlsx"), sheet = "ln_bhknots_second_revision"))
  
  # This is the covariance matrix provided by Linda from NJR in log_2_1
  rcs_second_names <- colnames(rcs_second_revision_mean)
  rcs_second_revision_covariance <- read_rcs_covariance(filename = paste0(data_directory,"/KNIPS Main input data.xlsx"),
                                                        sheetname = "second_revision_covariance",
                                                        par_names = rcs_second_names)
  
  rcs_second_revision <- mvrnorm(n_samples, mu = rcs_second_revision_mean,
                                 Sigma = rcs_second_revision_covariance)
  
  # Ensure it is a matrix
  if(n_samples == 1) {
    rcs_second_revision <- t(as.matrix(rcs_second_revision, nrow = 1))
  }
  
  
  # Transition 4
  # Post 1st revision to death
  # Same as transition 2
  
  # Transition 5
  # Post 2nd revision to death
  # Same as transition 2
  
  
  # Transition 6
  # Hazard rate of Post TKR to surgery death 
  # 90-day mortality from NJR
  mortality_90d_raw <- as.matrix(read_excel(paste0(data_directory,"/KNIPS Main input data.xlsx"), sheet = "90d_mortality"))
  mortality_row_index <- which(mortality_90d_raw[, "age_lower"] == age_range[1] & 
                                 mortality_90d_raw[, "age_upper"] == age_range[2] &
                                 mortality_90d_raw[, "gender"] == sample_gender)
  
  # Mortality probability is sampled from a beta distribution
  # Each sample is converted to a rate
  tkr_mortality_90d_beta_params <- as.numeric(mortality_90d_raw[mortality_row_index, "tkr_number"]) * c(
    as.numeric(mortality_90d_raw[mortality_row_index, "tkr_prop_dead"]),
    1 - as.numeric(mortality_90d_raw[mortality_row_index, "tkr_prop_dead"]))
  tkr_mortality_90d_probability <- rbeta(n_samples, shape1 = tkr_mortality_90d_beta_params[1], shape2 = tkr_mortality_90d_beta_params[2])
  tkr_mortality_90d <- as.matrix(-log(1 - tkr_mortality_90d_probability))
  
  # Rate is zero after 90 days
  tkr_mortality_90d <- cbind(tkr_mortality_90d, 0)
  colnames(tkr_mortality_90d) <- paste0("tkr_mortality_90d_", c(1, 2))
  tkr_mortality_times_90d <- c(0, 90/365.25)
  
  # Transition 7
  # Hazard rate of Post 1st revision to death
  # 90-day mortality from NJR
  
  # Mortality probability is sampled from a beta distribution
  # Each sample is converted to a rate
  revision_mortality_90d_beta_params <- as.numeric(mortality_90d_raw[mortality_row_index, "revision_number"]) * c(
    as.numeric(mortality_90d_raw[mortality_row_index, "revision_prop_dead"]),
    1 - as.numeric(mortality_90d_raw[mortality_row_index, "revision_prop_dead"]))
  revision_mortality_90d_probability <- rbeta(n_samples, shape1 = revision_mortality_90d_beta_params[1], shape2 = revision_mortality_90d_beta_params[2])
  revision_mortality_90d <- as.matrix(-log(1 - revision_mortality_90d_probability))
  
  # Rate is zero after 90 days
  revision_mortality_90d <- cbind(revision_mortality_90d, 0)
  colnames(revision_mortality_90d) <- paste0("revision_mortality_90d_", c(1, 2))
  revision_mortality_times_90d <- c(0, 90/365.25)
  
  # Transition 8
  # Hazard rate of Post 2nd revision to death
  rerevision_mortality_90d_beta_params <- as.numeric(mortality_90d_raw[mortality_row_index, "rerevision_number"]) * c(
    as.numeric(mortality_90d_raw[mortality_row_index, "rerevision_prop_dead"]),
    1 - as.numeric(mortality_90d_raw[mortality_row_index, "rerevision_prop_dead"]))
  rerevision_mortality_90d_probability <- rbeta(n_samples, shape1 = rerevision_mortality_90d_beta_params[1], shape2 = rerevision_mortality_90d_beta_params[2])
  rerevision_mortality_90d <- as.matrix(-log(1 - rerevision_mortality_90d_probability))
  
  # Beyond initial 90d period there is an ongoing risk of higher revision and surgical mortality
  # Probability of dying (within 90 days) times annual probability higher revision
  higher_revision_annual_mortality_probability <- rerevision_mortality_90d_probability * probability_higher_revision
  higher_revision_mortality <- as.matrix(-log(1 - higher_revision_annual_mortality_probability))
  
  # Initial 90d risk from 2nd revision, ongoing risk from further revision
  rerevision_mortality_90d <- cbind(rerevision_mortality_90d, higher_revision_mortality)
  colnames(rerevision_mortality_90d) <- paste0("rerevision_mortality_90d_", c(1, 2))
  rerevision_mortality_times_90d <- c(0, 90/365.25)
  
  ############################################################
  ## Costs and utilities
  ############################################################
  
  
  utilities_raw <- as.matrix(read_excel(paste0(data_directory,"/KNIPS Main input data.xlsx"), sheet = "utilities"))
  utility_row_index <- which(utilities_raw[, "age_lower"] == age_range[1] & 
                               utilities_raw[, "age_upper"] == age_range[2] &
                               utilities_raw[, "gender"] == sample_gender)
  
  
  revision_disutility_mean <- as.numeric(utilities_raw[utility_row_index, "disutilities"])
  revision_disutility_se <-  as.numeric(utilities_raw[utility_row_index, "disutilities SE"])
  utility_post_tkr_mean <-  as.numeric(utilities_raw[utility_row_index, "6 months after primary"])
  utility_post_1st_rev_mean <-  as.numeric(utilities_raw[utility_row_index, "6 months after revision"])
  utility_post_2nd_rev_mean <-   as.numeric(utilities_raw[utility_row_index, "6 months after revision"])
  utility_post_tkr_se <-  as.numeric(utilities_raw[utility_row_index, "6 months after primary SE"])
  utility_post_1st_rev_se <-  as.numeric(utilities_raw[utility_row_index, "6 months after revision SE"])
  utility_post_2nd_rev_se <-  as.numeric(utilities_raw[utility_row_index, "6 months after revision SE"])
  
  # Sample the utilities
  revision_disutility <- rnorm(n_samples, mean = revision_disutility_mean, sd = revision_disutility_se)
  utility_post_tkr <- rnorm(n_samples, mean = utility_post_tkr_mean, sd = utility_post_tkr_se)
  utility_post_1st_rev <- rnorm(n_samples, mean = utility_post_1st_rev_mean, sd = utility_post_1st_rev_se)
  utility_post_2nd_rev <- rnorm(n_samples, mean = utility_post_2nd_rev_mean, sd = utility_post_2nd_rev_se)
  
  # Implant costs 
  implant_costs_raw <- as.data.frame(read_excel(paste0(data_directory, "/KNIPS Main input data.xlsx"), sheet = "implant_costs"))
  
  implant_costs <- array(dim = c(n_strategies, n_samples), dimnames = list(implant_names, NULL))
  for(implant_name in implant_names) {
    implant_costs[implant_name, ] <- rnorm(n_samples, 
                                           mean = implant_costs_raw$mean[implant_costs_raw$implant_name == implant_name],
                                           sd = (implant_costs_raw$`95ci_ul`[implant_costs_raw$implant_name == implant_name] -
                                                   implant_costs_raw$`95ci_ll`[implant_costs_raw$implant_name == implant_name]) / (2 * 1.96))
  }
  
  # Load other costs
  other_costs_raw <- as.data.frame(read_excel(paste0(data_directory,"/KNIPS Main input data.xlsx"), sheet = "other_costs"))
  # Cost of primary TKR surgery
  tkr_surgery_cost <- other_costs_raw$value[other_costs_raw$parameter == "tkr_surgery_cost"]
  # Cost of revision surgery
  revision_cost <- other_costs_raw$value[other_costs_raw$parameter == "revision_cost"]
  
  ############################################################
  ## Return the inputs #######################################
  ############################################################
  
  # List of inputs for the economic model
  # This is a matrix of sample random parameters and some fixed parameters (e.g., knot locations, background mortality)
  model_inputs <- list()
  input_parameters <- matrix()
  
  # Fixed parameters
  model_inputs$ln_bhknots_first_revision <- ln_bhknots_first_revision
  model_inputs$ln_bhknots_second_revision <- ln_bhknots_second_revision
  model_inputs$mr <- mr
  model_inputs$tkr_mortality_times_90d <- tkr_mortality_times_90d
  model_inputs$revision_mortality_times_90d <- revision_mortality_times_90d
  model_inputs$rerevision_mortality_times_90d <- rerevision_mortality_times_90d
  
  model_inputs$tkr_surgery_cost <- tkr_surgery_cost
  model_inputs$revision_cost <- revision_cost
  
  
  # Random parameters
  input_parameters <- rcs_second_revision
  colnames(input_parameters) <- paste0("rcs_second_revision_", colnames(rcs_second_revision))
  input_parameters <- cbind(input_parameters, probability_higher_revision)
  for(implant_name in implant_names) {
    temp <- rcs_first_revision[[implant_name]]
    colnames(temp) <- paste0("rcs_first_revision_", implant_name, "_", colnames(rcs_first_revision[[implant_name]]))
    input_parameters <- cbind(input_parameters, temp)
  }
  
  input_parameters <- cbind(input_parameters, tkr_mortality_90d)
  input_parameters <- cbind(input_parameters, revision_mortality_90d)
  input_parameters <- cbind(input_parameters, rerevision_mortality_90d)
  
  input_parameters <- cbind(input_parameters, revision_disutility)
  input_parameters <- cbind(input_parameters, utility_post_tkr)
  input_parameters <- cbind(input_parameters, utility_post_1st_rev)
  input_parameters <- cbind(input_parameters, utility_post_2nd_rev)
  
  temp <- t(implant_costs)
  colnames(temp) <- paste0("implant_cost_", colnames(t(implant_costs)))
  input_parameters <- cbind(input_parameters, temp)

  # Add random parameters to the list
  model_inputs$inputs_parameters <- input_parameters
  
  return(model_inputs)
}