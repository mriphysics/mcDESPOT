#ifndef __signal_bSSFP180_B0_h__
#define __signal_bSSFP180_B0_h__

#include <Eigen/Core>
#include <limits>
#include "NaN.h"

class SSFP_steady_state_180_IV
{
  public:
    SSFP_steady_state_180_IV () :
      params (10)
    {
      params = Eigen::VectorXd::Constant(10,1,NaN);
    }

    void compute ();

    double& T1_F() { return params[0]; }
    double& T1_S() { return params[1]; }
    double& M0_F() { return params[2]; }
    double& M0_S() { return params[3]; }
    double& k_FS() { return params[4]; }
    double& k_SF() { return params[5]; }
    double& TR()   { return params[6]; }
    double& T2_F() { return params[7]; }
    double& T2_S() { return params[8]; }
    double& Delta() { return params[9]; }

    Eigen::VectorXd params;
    Eigen::VectorXd FA;

    Eigen::MatrixXd Mss;
    Eigen::VectorXd Mss_Sig;
};

#endif

