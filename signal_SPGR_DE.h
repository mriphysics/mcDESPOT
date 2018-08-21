#ifndef __signal_SPGR_DE_h__
#define __signal_SPGR_DE_h__

#include <Eigen/Core>
#include <limits>
#include "NaN.h"

class SPGR_steady_state_M0
{
  public:
    SPGR_steady_state_M0 () :
      params (5)
    {
      params = Eigen::VectorXd::Constant(5,1,NaN);
    }

    void compute ();

    double& T1_F() { return params[0]; }
    double& T1_S() { return params[1]; }
    double& M0_F() { return params[2]; }
    double& k_FS() { return params[3]; }
    double& TR() { return params[4]; }

    Eigen::VectorXd params;
    Eigen::VectorXd FA;

    Eigen::MatrixXd Mss;
    Eigen::VectorXd Mss_Sig;
};

#endif

