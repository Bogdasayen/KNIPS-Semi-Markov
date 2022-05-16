test1 <- rep(NA, 1200)
for(i in 1:1200) {
test1[i] <- Hsurvspline(
  x = 1,
  knots = ln_bhknots_second_revision,
  gamma = matrix(gamma_temp[i, ], nrow = 1),
  beta = matrix(beta_temp[i], nrow = 1, ncol = 1),
  X = yrs_to_1st_rev_temp$time_stop[i])
}
  

test2 <- Hsurvspline(
  x = 1, #matrix(1, nrow = 1200, ncol = 1),
  knots = ln_bhknots_second_revision,
  gamma = gamma_temp,
  beta = beta_temp * yrs_to_1st_rev_temp$time_stop,
  X = 1)#matrix(1:1200, ncol = 1))#matrix(yrs_to_1st_rev_temp$time_stop, ncol = 1))



Hsurvspline(
  x = c(1, 1, 1),
  knots = ln_bhknots_second_revision,
  gamma =  matrix(c(-2.522486, 0.8817235, -0.05348388, 0.08399305), nrow = 1),
  beta = matrix(-0.1367265, nrow = 1, ncol =1),
  X = c(1, 2, 3))

# Try to get first and second elements to match operation on 
# both simultaneously
Hsurvspline(
  x = c(1, 1),
  knots = ln_bhknots_second_revision,
  gamma = gamma_temp[c(1:2), ],
  beta = matrix(beta_temp[c(1:2)], ncol = 1), # changed to a row
  X = matrix(c(1, 2), nrow = 1))

# First element
Hsurvspline(
  x = c(1),
  knots = ln_bhknots_second_revision,
  gamma = matrix(gamma_temp[c(1), ], nrow = 1),
  beta = matrix(beta_temp[c(1)], nrow = 1, ncol = 1),
  X = matrix(c(1), nrow = 1, ncol =1))


