#ifndef __signal_bSSFP0_B0_DE_h__
#define __signal_bSSFP0_B0_DE_h__

#include <Eigen/Core>
#include <limits>
#include "NaN.h"

class SSFP_steady_state_IV
{
  public:
    SSFP_steady_state_IV () :
      params (8)
    {
      params = Eigen::VectorXd::Constant(8,1,NaN);
    }

    void compute ();

    double& T1_F() { return params[0]; }
    double& T1_S() { return params[1]; }
    double& M0_F() { return params[2]; }
    double& k_FS() { return params[3]; }
    double& TR()   { return params[4]; }
    double& T2_F() { return params[5]; }
    double& T2_S() { return params[6]; }
    double& Delta() { return params[7]; }

    Eigen::VectorXd params;
    Eigen::VectorXd FA;

    Eigen::MatrixXd Mss;
    Eigen::VectorXd Mss_Sig;
};

#endif

