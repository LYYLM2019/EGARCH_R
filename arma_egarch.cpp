#include <cmath>

#include "arma_egarch.h"

double var()
{
  double accum = 0.0;
  for_each(samples.begin(), samples.end(), [&](const double d) {
    accum += (d - mean) * (d - mean);
  });
  return accum / (data_size - 1);
}

EGARCH::EGARCH(const DATA &data):samples(data)
{
  data_size = samples.size();
  mean = accumulate(samples.begin(), samples.end(), 0) / data_size;
  variance = var();
  lower_bounds = -100*abs(m), -1, -1, -10, -10, -1, -10, 3;
  upper_bounds = 100*abs(m), 1, 1, 10, 10, 1, 10, 100;
  initial_parameters = abs(m), 0.1, 0.1, 1, 1, 0.1, 1, 10;
}

double EGARCH::egarchstdLLH(const PAR& parameters)
{
  double mu = parameters[0], ar = parameters[1], ma = parameters[2];
  double omega = parameters[3], alpha = parameters[4], bata = parameters[5];
  double gamma = parameters[6], shape = paramters[7];

  var = new DATA(data_size, 0);
  DATA log_var(data_size, 0), aux(data_size, 0), epsilon(data_size, 0), LLH(dta_size, 0);

  double fixvalue = log(tgamma(0.5 * (shape + 1)) / (tgamma(0.5 * shape) * (M_PI * (shape-2)) ^ 0.5));
  double expect = tgamma(0.5 * (shape - 1)) / tgamma(0.5 * shape) * (((shape-2) / M_PI) ^ 0.5);
  var[0] = variance;
  log_var[0] = log(variance);
  epsilon[0] = 0.001 * mean;
  z[0] = epsilon[0] / log_var[0];
  LLH[0] = -0.5 * log_var[0] - 0.5 * (1 + shape) * log(1 + (epsilon[0]^2 / (var[0] * (shape -2))));

  for(size_t i = 1; i < data_size; i++) {
    epsilon[i] = samples[i] - ar * samples[i-1] - ma * epsilon[i-1] - mu;
    log_var[i] = omega + alpha * z[i-1] + gamma * (abs(z[i-1]) - expect) + beta * log_var[i-1];
    var[i] = exp(log_var[i]);
    z[i] = epsilon[i] / (var[i]^0.5);
    LLH[i] = -0.5 * (log_var[i] - (1 + shape) * log(1 + (epsilon[i]^2 / (var[i] * (shape - 2)))))
  }
  objective = accumulate(LLH.begin(), LLH.end(), 0) + data_size * fixvalue;
  return objective;
}

Model EGARCH::fit()
{
  Model res;
  PAR final_paras = initial_parameters;
  find_min_box_constrained(bfgs_search_strategy(),
                           objective_delta_stop_strategy(1e-4),
                           egarchstdLLH, derivative(egarchstdLLH), final_paras, 0.1, 0.8)
  res.par = final_paras;
  res.simga2 = var;
  res.LLH = objective;
  return res;
}

EAGRCH::~EGARCH()
{
}

