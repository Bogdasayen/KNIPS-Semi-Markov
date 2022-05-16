# Second attempt at KNIPS model by modifying hesim_indivual_inhomogeneous_markov_1
# to match the inputs of knips_experiment_1.R

# V3 implements dependence on time to 1st revision
# V3a adds dependence of utilities and costs in 2nd revision state on time to 1st revision
# V3b optimises the survival calculations used for utilities and costs
# V4 Created a project and put input parameters in an Excel file

# Loop around line 469  calculating surival probabilities for 2nd revision dependent on time to first revision was a bottleneck
# Originally took 16 seconds (knips_experiment_3a.R) but reduced to 0.08 seconds (knips_experiment_3b.R)
# Overall this reduced runtime for 100pts/100samples from 4200 to 166 seconds
# For 1000pts/100samples took 2939 seconds

# Developed from hesim_indivual_inhomogeneous_markov_1.R
# hesim example of time inhomogeneous individual level Markov model
# Example is a hip replacement model
# https://hesim-dev.github.io/hesim/articles/markov-inhomogeneous-indiv.html

library(hesim)
library(data.table)
library(kableExtra)
library(flexsurv)
library(ggplot2)
library(MASS)
library(readxl)




# Load utility functions
source("code/knips_utils_1.R")

# Read in lifetables (from UK used in Coeliac screening - need to update)
lifetables <- read_excel("data/KNIPS Main input data.xlsx", sheet = "uk_lifetables")

############################################################
## Model specification
############################################################

n_samples <- 10
n_patients <- 10

discount_rate <- 0.035

# Not yet sure of best way to specify age - should it be a starting age and then match to a range?
age_range <- c(55,64)
starting_age <- ceiling(mean(age_range))
final_age <- 100

# Specify the gender
sample_gender <- "Female"

sys.time <-system.time({

# Implant details
implant_costs_raw <- as.data.frame(read_excel("data/KNIPS Main input data.xlsx", sheet = "implant_costs"))
n_implants <- dim(implant_costs_raw)[1]
implant_names <- implant_costs_raw$implant_name
  


# Time horizon is based on final age and starting age. The latter is age/gender dependent
time_horizon <- final_age - starting_age

# Number of time intervals only used for costs and utilities
# as these are discretized
n_time_intervals <- time_horizon


# Treatment strategies - one for each implant
strategies <- data.table(
  strategy_id = 1:n_implants,
  strategy_name = implant_names 
)
n_strategies <- nrow(strategies)

# Patients  
# Randomly sample their ages between the lower and upper age range limits

patients <- data.table(
  patient_id = 1:n_patients,
  gender = sample_gender,
  age = sample(c(age_range[1]:age_range[2]), n_patients, replace = TRUE)
)

# States - total hip replacement (THR) and revisions
states <- data.table( # Non-death health states
  state_id = 1:3,
  state_name = c("Post TKR", "Post 1st revision", "Post 2nd revision")
) 
n_states <- nrow(states)

# "hesim data"
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients, 
                        states = states)
print(hesim_dat)

# Nice labels - note that this includes death unless label set to NULL
labs <- get_labels(hesim_dat) #, death_label = NULL)
print(labs)

############################################################
## Transition rate parameters
############################################################


# Possible transitions
# Only 5 (first revision, 2nd revision, 3 death rates)
tmat <- rbind(c(NA, 1, NA, 2),
              c(NA, NA, 3, 4),
              c(NA, NA, NA, 5),
              c(NA, NA, NA, NA))
colnames(tmat) <- rownames(tmat) <- names(labs$state_id)
tmat

# Transition rates
# Transition 1
# From Post TKR to Post 1st revision

# Using all cases analysis from NJR
rcs_3knot_names <- c("cons", "rcs1", "rcs2", "rcs3")

implant_name <- "Cem CR_Fix Mono"

first_revision_filename <- paste0("data/", age_range[1], "-", age_range[2], "_", sample_gender, "_first_revision.xlsx")
rcs_first_revsision_mean_raw <- as.data.frame(read_excel(first_revision_filename, sheet = "rcs_first_revision_mean"))
ln_bhknots_first_revision_raw = as.data.frame(read_excel(first_revision_filename, sheet = "ln_bhknots_first_revision"))
rownames(rcs_first_revsision_mean_raw) <- rownames(ln_bhknots_first_revision_raw ) <- rcs_first_revsision_mean_raw[, 1]

rcs_first_revision_mean <- as.matrix(rcs_first_revsision_mean_raw[implant_name, -1])
ln_bhknots_first_revision <- as.matrix(ln_bhknots_first_revision_raw[implant_name, -1])

# This is the covariance matrix provided by Linda from NJR in log_2_1
rcs_first_revision_covariance <- read_rcs_covariance(filename = first_revision_filename,
                                                     sheetname = paste0(implant_name, "_cov"),
                                                 par_names = names(rcs_first_revision_mean))


rcs_first_revision <- mvrnorm(n_samples, mu = rcs_first_revision_mean,
                              Sigma = rcs_first_revision_covariance)

# Ensure it is a matrix
if(n_samples == 1) {
  rcs_first_revision <- t(as.matrix(rcs_first_revision, nrow = 1))
}


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
rcs_second_revision_mean <- as.matrix(read_excel("data/KNIPS Main input data.xlsx", sheet = "rcs_second_revision_mean"))
# Knots
ln_bhknots_second_revision = as.matrix(read_excel("data/KNIPS Main input data.xlsx", sheet = "ln_bhknots_second_revision"))

# This is the covariance matrix provided by Linda from NJR in log_2_1
rcs_second_names <- colnames(rcs_second_revision_mean)
rcs_second_revision_covariance <- read_rcs_covariance(filename = "data/KNIPS Main input data.xlsx",
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

############################################################
## Simulate disease progression
############################################################

# Create a list of disease progression results
disprog <- list()
# And of state occupancies
stateprobs <- list()

# Simulate one at a time
n_samples_temp <- 1

for(i_sample in 1:n_samples) {
  print(paste0("Sample ", i_sample))
  
  # Define the transition rate models
  
  # Simple function to convert vector to matrix with columns named constant ('cons')
  matrixv <- function(v, n = NULL){
    if (length(v) == 1) v <- rep(v, n_samples_temp) 
    m <- matrix(v)
    colnames(m) <- "cons"
    return(m)
  }
  
  # Function to generate random samples of underlying model parameters
  transmod_coef_def <- define_rng({
    
    ln_bhknots_first_revision <- fixed(ln_bhknots_first_revision[1,])

    ln_bhknots_second_revision = fixed(ln_bhknots_second_revision[1,])
    
    mr <- fixed(mr[,1])

    list(
      log_mr = lapply(as.list(log(mr)), matrixv),
      ln_bhknots_first_revision = as.matrix(ln_bhknots_first_revision),
      rcs_first_revision = rcs_first_revision[i_sample, ],
      rcs_second_revision = rcs_second_revision[i_sample, ],
      ln_bhknots_second_revision = as.matrix(ln_bhknots_second_revision)
    )
  }, n = n_samples_temp)
  transmod_coef <- eval_rng(transmod_coef_def)
  
  # Create a list of randomly sampled survival models for transitions
  # This uses parameter samples in transmod_coef and specified distributions
  # First ensure the coefficients are matrices
  if(n_samples_temp == 1) {
    transmod_coef$rcs_first_revision <- t(as.matrix(transmod_coef$rcs_first_revision, nrow = 1))
    transmod_coef$rcs_second_revision <- t(as.matrix(transmod_coef$rcs_second_revision, nrow = 1))
  }
  transmod_params <- params_surv_list(
    # 1. Post TKR to Post 1st revision
    params_surv(coefs = list(gamma = transmod_coef$rcs_first_revision),
                aux = list(knots = transmod_coef$ln_bhknots_first_revision,
                           scale = "log_hazard",
                           timescale = "log"),
                dist = "survspline"),
    
    
    # 2. Post TKR to death
    params_surv(coefs = transmod_coef$log_mr,
                aux = list(time = c(0:time_horizon)),
                dist = "pwexp"),
    
    # 3. Post 1st revision to post 2nd revision
    # Includes dependence on time from TKR to 1st revision
    params_surv(coefs = list(gamma = transmod_coef$rcs_second_revision),
                aux = list(knots = transmod_coef$ln_bhknots_second_revision,
                           scale = "log_hazard",
                           timescale = "log"),
                dist = "survspline"),
    
    # 4. Post 1st revision to death
    params_surv(coefs = transmod_coef$log_mr,
                aux = list(time = c(0:time_horizon)),
                dist = "pwexp"),
    
    # 5. Post 2nd revision to death
    params_surv(coefs = transmod_coef$log_mr,
                aux = list(time = c(0:time_horizon)),
                dist = "pwexp")
  )
  
  # Now simulate the outcomes for all patients for this sample
  
  transmod_data <- expand(hesim_dat, by = c("strategies", "patients"))
  head(transmod_data)
  # Add numerical covariates
  # Need to add constants to match the spline parameters
  transmod_data[, rcs1 := 1] 
  transmod_data[, rcs2 := 1]
  transmod_data[, rcs3 := 1]
  transmod_data[, cons := 1]
  transmod_data[, yrs_to_1st_rev := 0]
  transmod_data[, male := ifelse(gender == "Male", 1, 0)]

  
  # Transition model
  transmod <- create_IndivCtstmTrans(transmod_params, 
                                     input_data = transmod_data,
                                     trans_mat = tmat,
                                     clock = "mix",
                                     reset_states = c(1, 2, 3),
                                     start_age = patients$age)
  
  
  set.seed(2243534)
  disprog[[i_sample]] <- transmod$sim_disease(max_t = starting_age, max_age = (starting_age + time_horizon))
  
  # Update the data using results of this initial simulation
  # Use the time to transition from state 1 to state 2
  # Not all patients on all transitions actually move from 1 to 2
  disprog_temp <- disprog[[i_sample]][from == 1 & to == 2, ]
  
  # Using slow loop but might be better way
  for(i_transition in 1:dim(disprog_temp)[1]) {
    transmod_data[transmod_data$patient_id == as.numeric(disprog_temp[i_transition, "patient_id"]) &
                    transmod_data$strategy_id == as.numeric(disprog_temp[i_transition, "strategy_id"]), "yrs_to_1st_rev"] <-
      disprog_temp[i_transition, "time_stop"]
    
  }
  
  # Re-run simulation with the (correct) time to first revision
  # Transition model
  transmod <- create_IndivCtstmTrans(transmod_params, 
                                     input_data = transmod_data,
                                     trans_mat = tmat,
                                     clock = "mix",
                                     reset_states = c(1, 2, 3),
                                     start_age = patients$age)
  set.seed(2243534)
  disprog[[i_sample]] <- transmod$sim_disease(max_t = starting_age, max_age = (starting_age + time_horizon))
  # Set the sample number
  disprog[[i_sample]]$sample <- i_sample
  
  # Also simulate average state probabilities
  stateprobs[[i_sample]] <- transmod$sim_stateprobs(t = 0:(time_horizon))
  stateprobs[[i_sample]]$sample <- i_sample
  
  
} # End loop over samples


# Create disprog_combined
# I'd recommend a slight tweak in that I'd use the IndivCtstmTrans class to 
# simulate disease via $sim_disease(). I'd store each resulting "disprog" object 
# in a list and then row bind them together (e.g., via rbindlist()). Once you've 
# binded them you need to make sure two things are correct: (i) the attributes 
# and (ii) the sorting order. Use setattr() to make sure the size and absorbing 
# attributes are defined correctly. The table should be sorted by sample, 
# strategy_id, patient_id, and time_start (e.g., via setorderv()).

disprog_combined <-rbindlist(disprog)
stateprobs_combined <- rbindlist(stateprobs)

# Correctly set the size attributes
# Check if using attributes(disprog_combined)
attributes_size <- c(n_samples, n_strategies, n_patients, n_states + 1)
names(attributes_size) <- c("n_samples", "n_strategies", "n_patients", "n_states")
setattr(disprog_combined, "size", attributes_size)



############################################################
## Costs and utilities
############################################################


utilities_raw <- as.matrix(read_excel("data/KNIPS Main input data.xlsx", sheet = "utilities"))
rownames(utilities_raw) <- utilities_raw[, "parameter"]

revision_disutility_mean <- as.numeric(utilities_raw["revision_disutility_mean", "value"])
revision_disutility_se <-  as.numeric(utilities_raw["revision_disutility_se", "value"])
utility_post_tkr_mean <-  as.numeric(utilities_raw["utility_post_tkr_mean", "value"])
utility_post_1st_rev_mean <-  as.numeric(utilities_raw["utility_post_1st_rev_mean", "value"])
utility_post_2nd_rev_mean <-  as.numeric(utilities_raw["utility_post_2nd_rev_mean", "value"])
utility_post_tkr_se <-  as.numeric(utilities_raw["utility_post_tkr_se", "value"])
utility_post_1st_rev_se <-  as.numeric(utilities_raw["utility_post_1st_rev_se", "value"])
utility_post_2nd_rev_se <-  as.numeric(utilities_raw["utility_post_2nd_rev_se", "value"])

# Sample the utilities
revision_disutility <- rnorm(n_samples, mean = revision_disutility_mean, sd = revision_disutility_se)
utility_post_tkr <- rnorm(n_samples, mean = utility_post_tkr_mean, sd = utility_post_tkr_se)
utility_post_1st_rev <- rnorm(n_samples, mean = utility_post_1st_rev_mean, sd = utility_post_1st_rev_se)
utility_post_2nd_rev <- rnorm(n_samples, mean = utility_post_2nd_rev_mean, sd = utility_post_2nd_rev_se)

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

# Intervals are 1 year if time_horizon is equal to n_time_intervals
for(time_start_ in 0:(n_time_intervals - 1)) {
  # Currently the same for all strategies but will eventually depend on them
  # by changing the rcs_first_revision and ln_bhknots_first_revision
  probability_1st_revision[probability_1st_revision$time_start == time_start_, "value"] <- 
    rep(exp(-Hsurvspline(
      x = time_start_,
      gamma = rcs_first_revision,
      knots = ln_bhknots_first_revision)) - 
        exp(-Hsurvspline(
          x = time_start_ + time_horizon/n_time_intervals,
          gamma = rcs_first_revision,
          knots = ln_bhknots_first_revision
        )), n_strategies)
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
for(time_start_ in 0:(n_time_intervals - 1)) {
  print(paste0("Time interval ", time_start_, "/", n_time_intervals))
  
  gamma_temp <- rcs_second_revision[, -2]
  beta_temp <- rcs_second_revision[, 2]
  
  # Create duplicates for each patient and strategy
  # These only vary by sample
  gamma_temp <- gamma_temp[rep(c(1:n_samples), n_patients * n_strategies), ]
  beta_temp <- rep(beta_temp, n_patients * n_strategies)
  
  
  probability_2nd_revision[probability_2nd_revision$time_start == time_start_, "value"] <- 
    exp(-Hsurvspline(
      x = time_start_,
      knots = ln_bhknots_second_revision,
      gamma = gamma_temp,
      beta = beta_temp,
      X = yrs_to_1st_rev_temp$time_stop)) - 
    exp(-Hsurvspline(
      x = time_start_ + time_horizon/n_time_intervals,
      knots = ln_bhknots_second_revision,
      gamma = gamma_temp,
      beta = beta_temp, 
      X = yrs_to_1st_rev_temp$time_stop)) 
  
} # End loop over times





# Utilities depend on the probabilities of revision
# This may not work if n_patients or n_samples becomes too big
utility_tbl <- stateval_tbl(
  data.table(strategy_id = rep(strategies$strategy_id, each = (n_patients * n_states * n_samples * n_time_intervals)),
             patient_id = rep(rep(1:n_patients, times = n_strategies), each = (n_states * n_samples * n_time_intervals)),
             state_id = rep(rep(states$state_id, times = n_strategies * n_patients), each = (n_samples * n_time_intervals)),
             sample = rep(rep(1:n_samples, times = n_strategies * n_patients * n_states), each = n_time_intervals),
             time_start = rep(0:(n_time_intervals - 1), times = (n_strategies * n_patients * n_states * n_samples)),
             value = rep(c(0, 0, 0), times = n_strategies * n_patients * n_samples * n_time_intervals)),
  dist = "custom"
)


# Utilities in Post TKR are dependent on time (clock-forward as all start in Post TKR)
for(time_start_ in 0:(n_time_intervals - 1)) {
  # Currently the same for all strategies but will eventually depend on them
  # by changing the rcs_first_revision and ln_bhknots_first_revision
  utility_tbl[utility_tbl$state_id == 1 & utility_tbl$time_start == time_start_, "value"] <- rep(utility_post_tkr, n_strategies * n_patients) -# general utility
    # Disutility times probability of revision during this time interval
    + rep(revision_disutility, n_strategies * n_patients) *
    rep(unlist(probability_1st_revision[probability_1st_revision$time_start == time_start_, "value"]), times = n_patients)
}




# Utilities in Post 1st revision are dependent on time 
# Note that this is clock reset
for(time_start_ in 0:(n_time_intervals - 1)) {
  # Currently the same for all strategies but will eventually depend on them
  # by changing the rcs_first_revision and ln_bhknots_first_revision
  utility_tbl[utility_tbl$state_id == 2 & utility_tbl$time_start == time_start_, "value"] <- rep(utility_post_1st_rev, n_strategies * n_patients) -# general utility
    # Disutility times probability of revision during this time intervals
    + rep(revision_disutility, n_strategies * n_patients) * 
    probability_2nd_revision[probability_2nd_revision$time_start == time_start_, "value"]
  
}

# Utility post 2nd revision the same for all strategies
# Need to change so time dependent
# and dependent on probabilities of 3rd or higher revision
utility_tbl[utility_tbl$state_id == 3, "value"] <- rep(utility_post_2nd_rev, each = n_strategies * n_patients * n_time_intervals)


# Check that all strategies currently have the same utility in Post TKR
#utility_tbl[utility_tbl$state_id == 1 & utility_tbl$time_start == 20 &
#             utility_tbl$sample == 1, ]
# And in Post 1st revision
#utility_tbl[utility_tbl$state_id == 2 & utility_tbl$time_start == 20 &
#     utility_tbl$sample == 1, ]


head(utility_tbl)

# Implant costs 
implant_costs <- array(dim = c(n_strategies, n_samples), dimnames = list(implant_names, NULL))
for(implant_name in implant_names) {
  implant_costs[implant_name, ] <- rnorm(n_samples, 
                                         mean = implant_costs_raw$mean[implant_costs_raw$implant_name == implant_name],
                                         sd = (implant_costs_raw$`95ci_ul`[implant_costs_raw$implant_name == implant_name] -
                                                 implant_costs_raw$`95ci_ll`[implant_costs_raw$implant_name == implant_name]) / (2 * 1.96))
}

# Load other costs
other_costs_raw <- as.data.frame(read_excel("data/KNIPS Main input data.xlsx", sheet = "other_costs"))
# Cost of primary TKR surgery
tkr_surgery_cost <- other_costs_raw$value[other_costs_raw$parameter == "tkr_surgery_cost"]
# Cost of revision surgery
revision_cost <- other_costs_raw$value[other_costs_raw$parameter == "revision_cost"]



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
             value = rep(c(0, 0, 0), times = n_strategies * n_patients * n_samples * n_time_intervals)),
  dist = "custom"
)

# Costs in Post TKR are dependent on time (clock-forward as all start in Post TKR)
for(time_start_ in 0:(n_time_intervals - 1)) {
  # Currently the same for all strategies but will eventually depend on them
  # by changing the rcs_first_revision and ln_bhknots_first_revision
  medcost_tbl[medcost_tbl$state_id == 1 & medcost_tbl$time_start == time_start_, "value"] <-  
    # Surgery cost times probability of revision during this time intervals
    + rep(revision_cost, n_strategies * n_patients) *
  rep(unlist(probability_1st_revision[probability_1st_revision$time_start == time_start_, "value"]), times = n_patients)
  
}

# Costs in Post 1st revision are dependent on time 
# Note this is clock reset
for(time_start_ in 0:(n_time_intervals - 1)) {
  # Currently the same for all strategies but will eventually depend on them
  # by changing the rcs_first_revision and ln_bhknots_first_revision
  medcost_tbl[medcost_tbl$state_id == 2 & medcost_tbl$time_start == time_start_, "value"] <- 
    # Cost times probability of revision during this time intervals
    + rep(revision_cost, n_strategies * n_patients) * 
    probability_2nd_revision[probability_2nd_revision$time_start == time_start_, "value"]
}

# Cost post 2nd revision the same for all strategies
# Need to change so time dependent
# and dependent on probabilities of 3rd or higher revision
medcost_tbl[medcost_tbl$state_id == 3, "value"] <- rep(revision_cost, each = n_samples * n_strategies * n_patients * n_time_intervals) * 
  # Using mean probability of 2nd revision for now but this needs to be updated
  mean(unlist(lapply(probability_2nd_revision[, "value"], as.numeric)))


############################################################
## Simulate costs and QALYs
############################################################

# Utility and cost models
# Utility
utilitymod <- create_StateVals(utility_tbl, n = n_samples, 
                               hesim_data = hesim_dat)

# Costs
# The 'starting' option means costs are only applied when patients enter the state
implantcostmod <- create_StateVals(implantcost_tbl, n = n_samples,
                                   method = "starting", hesim_data = hesim_dat)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples,
                               hesim_data = hesim_dat)
costmods <- list(Drug = implantcostmod,
                 Medical = medcostmod)

# Economic model combining everything
econmod <- IndivCtstm$new(trans_model = transmod,
                          utility_model = utilitymod,
                          cost_models = costmods)

# Use the combined disease progression results
econmod$disprog_ <- disprog_combined

# Look at simulations to ensure the combination was correct
head(econmod$disprog_)

# And the combined state occupancy probabilities
econmod$stateprobs_ <- stateprobs_combined

econmod$stateprobs_
# Average probability of occupying Post TKR
with(econmod$stateprobs_, mean(prob[state_id == 1]))
# Average in post 1st revision
with(econmod$stateprobs_, mean(prob[state_id == 2]))
# Average in post 2nd revision
with(econmod$stateprobs_, mean(prob[state_id == 3]))
# Average probability of occupying death
with(econmod$stateprobs_, mean(prob[state_id == 4]))

# Simulate QALYs and costs with discount rates
# 3.5% discount rate
econmod$sim_qalys(dr = discount_rate)
econmod$sim_costs(dr = discount_rate)


############################################################
## Analyse results
############################################################


ce_sim <- econmod$summarize()
summary(ce_sim, labels = labs) %>%
  format()

# And calculate ICER
cea_pw_out <- cea_pw(ce_sim, comparator = 1, dr_qalys = discount_rate, dr_costs = discount_rate,
                     k = seq(0, 25000, 500))
icer(cea_pw_out, labels = labs) %>%
  format(digits_qalys = 3)

}) # end system.time