#ifndef __signal_SPGR_SP_h__
#define __signal_SPGR_SP_h__

#include <Eigen/Core>
#include <limits>
#include "NaN.h"

class SPGR_SP_SS
{
  public:
    SPGR_SP_SS () :
      params (3)
    {
      params = Eigen::VectorXd::Constant(3,1,NaN);
    }

    void compute ();

    double& T1() { return params[0]; }
    double& M0() { return params[1]; }
    double& TR() { return params[2]; }

    Eigen::VectorXd params;
    Eigen::VectorXd FA;

    //Eigen::VectorXd Mss;
    Eigen::VectorXd Mss_Sig;
};

#endif
