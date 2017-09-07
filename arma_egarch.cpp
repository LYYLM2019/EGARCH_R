#include <cmath>
#include <numeric>
#include <fstream>
#include <string>
#include <limits>
#include <cfloat>

#include "arma_egarch.h"

typedef matrix<double,0,1> column_vector;
double rosen (const column_vector& m)
{
  const double x = m(0);
  const double y = m(1);
  return 100.0*pow(y - x*x,2) + pow(1 - x,2);
}

/*double EGARCH::var_calculation()
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
  variance = var_calculation();
  lower_bounds.set_size(8);
  upper_bounds.set_size(8);
  initial_parameters.set_size(8);
  lower_bounds = -100*abs(mean), -1.0, -1.0, -10.0, -10.0, -1.0, -10.0, 3.0;
  upper_bounds = 100*abs(mean), 1.0, 1.0, 10.0, 10.0, 1.0, 10.0, 100.0;
  initial_parameters = abs(mean), 0.1, 0.1, 1.0, 1.0, 0.1, 1.0, 100.0;
}

double EGARCH::egarchstdLLH(const PAR& parameters)
{
  double mu = parameters(0), ar = parameters(1), ma = parameters(2);
  double omega = parameters(3), alpha = parameters(4), beta = parameters(5);
  double gamma = parameters(6), shape = parameters(7);

  var = DATA(data_size, 0);
  DATA log_var(data_size, 0), aux(data_size, 0), epsilon(data_size, 0), LLH(data_size, 0);

  double fixvalue = log(tgamma(0.5 * (shape + 1)) / (tgamma(0.5 * shape) * sqrt(M_PI * (shape-2))));
  double expect = tgamma(0.5 * (shape - 1)) / tgamma(0.5 * shape) * (sqrt((shape-2) / M_PI));
  var[0] = variance;
  log_var[0] = log(variance);
  epsilon[0] = 0.001 * mean;
  aux[0] = epsilon[0] / log_var[0];
  LLH[0] = -0.5 * log_var[0] - 0.5 * (1 + shape) * log(1 + (epsilon[0] * epsilon[0] / (var[0] * (shape -2))));

  for(size_t i = 1; i < data_size; i++) {
    epsilon[i] = samples[i] - ar * samples[i-1] - ma * epsilon[i-1] - mu;
    log_var[i] = omega + alpha * aux[i-1] + gamma * (abs(aux[i-1]) - expect) + beta * log_var[i-1];
    var[i] = exp(log_var[i]);
    aux[i] = epsilon[i] / sqrt(var[i]);
    LLH[i] = -0.5 * (log_var[i] - (1 + shape) * log(1 + (epsilon[i] * epsilon[i] / (var[i] * (shape - 2)))));
  }
  objective = accumulate(LLH.begin(), LLH.end(), 0) + data_size * fixvalue;
  return objective;
}

Model EGARCH::fit()
{
  Model res;
  //PAR final_paras = initial_parameters;
  column_vector final_paras(2);
  final_paras = 0.1, 0.1;
  find_min_box_constrained(bfgs_search_strategy(),
                           objective_delta_stop_strategy(1e-4),
                           rosen, derivative(rosen), final_paras, 0.1, 0.8);
  //res.par = final_paras;
  //res.sigma2 = var;
  //res.LLH = objective;
  cout << final_paras;
  return res;
}

EGARCH::~EGARCH()
{
}
*/
void fit()
{
  column_vector final_paras(2);
  final_paras = 0.1, 0.1;
  find_min_box_constrained(bfgs_search_strategy(),
                           objective_delta_stop_strategy(1e-4),
                           rosen, derivative(rosen), final_paras, 0.1, 0.8);
  cout << final_paras;
}

void read_data(DATA& dt)
{
  ifstream data("d4.txt");
  string line;
  while(getline(data, line)) {
    dt.push_back(stod(line));
  }
}

unsigned int data_size;
DATA samples;
double variance1;
double mean1;
double objective;

double var_calculation()
{
  double accum = 0.0;
  for_each(samples.begin(), samples.end(), [&](const double d) {
    accum += (d - mean1) * (d - mean1);
  });
  return accum / (data_size - 1);
}

double fRand(double fMin, double fMax)
{
  double f = (double)std::rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

double egarchstdLLH(const PAR& parameters)
{
  double mu = parameters(0), ar = parameters(1), ma = parameters(2);
  double omega = parameters(3), alpha = parameters(4), beta = parameters(5);
  double gamma = parameters(6), shape = parameters(7);

  DATA var = DATA(data_size, 0.0);
  DATA log_var(data_size, 0.0), aux(data_size, 0.0), epsilon(data_size, 0.0), LLH(data_size, 0.0);

  double fixvalue = log(tgamma(0.5 * (shape + 1)) / (tgamma(0.5 * shape) * sqrt(M_PI * (shape - 2))));
  double expect = tgamma(0.5 * (shape - 1)) / tgamma(0.5 * shape) * (sqrt((shape - 2) / M_PI));
  var[0] = variance1;
  log_var[0] = log(variance1);
  epsilon[0] = 0.001 * mean1;
  aux[0] = epsilon[0] / sqrt(var[0]);
  LLH[0] = fixvalue - 0.5 * log_var[0] - 0.5 * (1 + shape) * log(1 + (epsilon[0] * epsilon[0] / (var[0] * (shape -2))));
  //cout << LLH[0] << endl;
  for(size_t i = 1; i < data_size; i++) {
    epsilon[i] = samples[i] - ar * samples[i-1] - ma * epsilon[i-1] - mu;
    log_var[i] = omega + alpha * aux[i-1] + gamma * (abs(aux[i-1]) - expect) + beta * log_var[i-1];
    if (isnan(log_var[i]) || log_var[i] < -1000000) {
      cout << parameters << endl;
      return DBL_MAX;
    }
    var[i] = exp(log_var[i]);
    if (isnan(var[i]) || var[i] + 1.0 == var[i]) {
 //     cout << parameters << endl;
      return 1e+20;
      //return std::numeric_limits<double>::infinity();
    }
    aux[i] = epsilon[i] / sqrt(var[i]);
    if (isnan(aux[i]) || aux[i] + 1.0 == aux[i])
      return 1e+10;
    LLH[i] = fixvalue - 0.5 * log_var[i] - 0.5 * (1 + shape) * log(1 + (epsilon[i] * epsilon[i] / (var[i] * (shape - 2))));
    //cout << log_var[i] << ' ' << var[i] << ' ' << aux[i] << ' ' << LLH[i] << endl;// ' ' << LLH[i] << endl;
  }
  objective = accumulate(LLH.begin(), LLH.end(), 0.0);// + data_size * fixvalue;
  cout << "obj: " << objective << endl;
  return -objective;
}

int main() {
  read_data(samples);
  data_size = samples.size();
  cout << "sample size: " << samples.size() << endl;
  cout << "sample sum: " << accumulate(samples.begin(), samples.end(), 0.0) << endl;
  mean1 = accumulate(samples.begin(), samples.end(), 0.0) / data_size;
  cout << "sample mean: " << mean1 << endl;
  variance1 = var_calculation();
  cout << "sample variance: " << variance1 << endl;
  PAR lower_bounds;
  PAR upper_bounds;
  PAR initial_parameters;
  lower_bounds.set_size(8);
  upper_bounds.set_size(8);
  initial_parameters.set_size(8);
  lower_bounds = -100*abs(mean1), -1.0, -1.0, -10.0, -10.0, -1.0, -10.0, 3;
  upper_bounds = 100*abs(mean1), 1.0, 1.0, 10.0, 10.0, 1.0, 10.0, 100.0;
  initial_parameters = abs(mean1), 0.02, 0.02, 0.01, 0.02, 0.3, 0.03, 5.0;
  //cout << "initial parameters: " << initial_parameters << endl;
  //cout << "initial LLH: " << egarchstdLLH(initial_parameters) << endl;
  //cout << "derivative result: " << derivative(egarchstdLLH)(initial_parameters) << endl;
  /*cout << find_max(bfgs_search_strategy(),
                   objective_delta_stop_strategy(0.1).be_verbose(),
                   egarchstdLLH, derivative(egarchstdLLH),
                   initial_parameters, 2500.0) << endl;
  */
  cout << find_min_box_constrained(bfgs_search_strategy(),
                                   objective_delta_stop_strategy(1e-8).be_verbose(),
                                   egarchstdLLH, derivative(egarchstdLLH),
                                   initial_parameters, lower_bounds, upper_bounds) << endl;
  //find_max_using_approximate_derivatives(bfgs_search_strategy(), objective_delta_stop_strategy(1e-6), egarchstdLLH, initial_parameters, 500.0);
  cout << "result: \n" << initial_parameters << endl;
  PAR s;
  s.set_size(8);
  //s = 0.0017785, 0.0199709, 0.0199651, -7.10829, 0.02025, 0.570112, 5.25887, 4.43677;
  //s = 0.00169982, 1, 1, 10, 10, 1, 10, 100;
  s = -2.172466e-05,-3.183580e-02,-1.666200e-01,-2.277217e+00,-6.394999e-02,8.436640e-01,3.188545e-01,5.671221e+00;
  cout << "LLH:," << egarchstdLLH(s) << endl;
  return 0;
}
