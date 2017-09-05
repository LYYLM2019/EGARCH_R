# global variable
samples <- numeric(100)
len <- length(samples)

egarchstdLLH <- function (parameters) {
  fixvalue <- log(gamma(0.5 * (parameters[8] + 1)) / (gamma(0.5 * parameters[8]) * (pi * (parameters[8]-2)) ^ 0.5))
  expect = gamma(0.5 * (parameters[8] - 1)) / gamma(0.5 * parameters[8]) * (((parameters[8]-2) / pi) ^ 0.5)
  sigma2<- numeric(len)
  log_sigma2 <- numeric(len)
  epsilon <- numeric(len)
  LLH <- numeric(len)
  z <- numeric(len)

  sigma2[1] <- var(samples)
  log_sigma2[1] <- log(sigma2[1])
  epsilon[1] <- 0.001 * mean(samples)
  z[1] <- epsilon[1] / (sigma2[1]^0.5)
  LLH[1]<- fixvalue - 0.5 * log_sigma2[1] - 0.5 * (1 + parameters[8]) * log(1 + (epsilon[1]^2 / (sigma2[1] * (parameters[8] -2))))

  for (i in 2:len) {
    epsilon[i] = samples[i] - parameters[2] * samples[i-1] - parameters[3] * epsilon[i-1] - parameters[1]
    log_sigma2[i] = parameters[4] + parameters[5] * z[i-1] + parameters[7] * (abs(z[i-1]) - expect) + parameters[6] * log_sigma2[i-1]
    sigma2[i] = exp(log_sigma2[i])
    z[i] = epsilon[i] / (sigma2[i]^0.5)
    LLH[i] <- fixvalue - 0.5 * log_sigma2[i] - 0.5 * (1 + parameters[8]) * log(1 + (epsilon[i]^2 / (sigma2[i] * (parameters[8] - 2))))
    print(LLH[i])
  }
  return (-sum(LLH))
}


egarchfit <- function (samples, method = 'solnp') {
  m=mean(samples)
  v=var(samples)
  parameters <- c(m, 0.02, 0.02, 0.01, 0.02, 0.3, 0.03, 5)
  lowerBounds <- c(100*m, -1, -1, -10, -10, -1, -10, 3)
  upperBounds <- c(-100*m, 1, 1, 10, 10, 1, 10, 100)
  if (method == 'solnp')
    model = solnp(parameters, egarchstdLLH, LB=lowerBounds, UB=upperBounds)
  else if (method == 'nlminb')
    nlminb(parameters, egarchstdLLH, lower=lowerBounds, upper=upperBounds)
}
