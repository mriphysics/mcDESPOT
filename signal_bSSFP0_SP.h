#ifndef __signal_bSSFP0_SP_h__
#define __signal_bSSFP0_SP_h__

#include <Eigen/Core>
#include <limits>
#include "NaN.h"

class SSFP_SP_SS
{
  public:
    SSFP_SP_SS () :
      params (5)
    {
      params = Eigen::VectorXd::Constant(5,1,NaN);
    }

    void compute ();

    double& T1() { return params[0]; }
    double& M0() { return params[1]; }
    double& TR() { return params[2]; }
    double& T2() { return params[3]; }
    double& Delta() { return params[4]; }

    Eigen::VectorXd params;
    Eigen::VectorXd FA;

    Eigen::MatrixXd Mss;
    Eigen::VectorXd Mss_Sig;
};

#endif
