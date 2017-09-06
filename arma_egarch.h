/***************************************************************
 * Copyright (C) COINUT PTE. LTD. - All Rights Reserved
 *
 * Unauthorized copying of this file, via any medium is strictly
 * prohibited Proprietary and confidential.
 *
 * Written by Jingli Cai <jingli@coinut.com>, 2017 - 2018
 *
 * ARMA(1,1)-EGARCH(1,1) MODEL
 *
 ***************************************************************/

#ifndef ARMA_EGARCH_H
#define ARMA_EGARCH_H

#include <iostream>
#include <vector>
#include <dlib>

using namespace std;
using namespace dlib;

typedef matrix<double, 0, 8> PAR;
typedef vector<double> DATA;

struct Model {
  PAR par;
  double LLH;
  DATA sigma2;
};

class EGARCH {
  public:
    EGARCH(const DATA &data):samples(data) {
    }
    double egarchstdLLH(const PAR& parameters);
    Model fit();
    ~EGARCH();
  private:
    DATA samples;
    DATA var;
    double objective;
    unsigned int data_size;
    double mean, variance;
    PAR uppper_bounds, lower_bounds;
    PAR initial_parameters;
    Model model;
};

#endif
