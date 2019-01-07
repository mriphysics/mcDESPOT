#ifndef __signal_SPGR_h__
#define __signal_SPGR_h__

#include <Eigen/Core>
#include <limits>
#include "NaN.h"

class SPGR_steady_state_M0
{
  public:
    SPGR_steady_state_M0 () :
      params (7)
    {
      params = Eigen::VectorXd::Constant(7,1,NaN);
    }

    void compute ();

    double& T1_F() { return params[0]; }
    double& T1_S() { return params[1]; }
    double& M0_F() { return params[2]; }
    double& M0_S() { return params[3]; }
    double& k_FS() { return params[4]; }
    double& k_SF() { return params[5]; }
    double& TR() { return params[6]; }

    Eigen::VectorXd params;
    Eigen::VectorXd FA;

    Eigen::MatrixXd Mss;
    Eigen::VectorXd Mss_Sig;
};

#endif

