#include <unsupported/Eigen/MatrixFunctions>
#include "signal_SPGR_SP.h"
#include <iostream>
#include <cmath>

using namespace Eigen;

void SPGR_SP_SS::compute ()
{

  // Calculate remaining parameters.
  double R1 = 1.0/T1();

  // Define Bloch-McConnell terms:
  double C = (R1 * M0());

  double A = -R1 ;

  Mss_Sig.resize (FA.size());

  double AinvC = C/A;

  A *= TR();
  double em = exp(A);

  for (int n = 0; n < FA.size(); n++) {

    double T = std::cos(FA(n));

    Mss_Sig(n) = std::sin(FA(n)) * pow((1 - (em * T)),-1) * ((em - 1) * AinvC);

  }

}
