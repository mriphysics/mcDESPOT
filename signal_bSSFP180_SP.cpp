#include <unsupported/Eigen/MatrixFunctions>
#include "signal_bSSFP180_SP.h"
#include <complex>
#include <iostream>

using namespace Eigen;

void SSFP_SP_SS_180::compute ()
{
  // Calculate remaining parameters.
  double R1 = 1.0/T1();
  double R2 = 1.0/T2();

  // Define Bloch-McConnell terms:
  VectorXd C(3);
  C << 0.0, 0.0, (R1 * M0());

  MatrixXd A(3,3);
  A << -R2, 0.0, 0.0,
       0.0, -R2, 0.0, 
       0.0, 0.0, -R1;
    
  Mss.resize (FA.size(), 3);
  Mss_Sig.resize (FA.size());

  PartialPivLU<MatrixXd> lu (A);
  VectorXd AinvC(3);
  AinvC = lu.solve (C);

  A *= TR();
  MatrixXd em(3,3);
  em = A.exp();

  const std::complex<double> i(0.0,1.0); 

  double PCTheta = 3.14159265358979323846 + Delta();

  MatrixXd PCM(3,3);

  PCM << std::cos(PCTheta), -std::sin(PCTheta), 0.0,
         std::sin(PCTheta), std::cos(PCTheta), 0.0, 
         0.0, 0.0, 1.0;

  MatrixXd T(3,3);

  for (int n = 0; n < FA.size(); n++) {

    T << 1.0, 0.0, 0.0,
         0.0, std::cos(FA(n)), std::sin(FA(n)),
         0.0, -std::sin(FA(n)), std::cos(FA(n));

    lu.compute (PCM - (em * T));
    Mss.row(n) = lu.solve ((em - MatrixXd::Identity(3,3)) * AinvC);

    // Extract signal component.
    Mss_Sig(n) = std::abs((Mss(n,0) + (i * Mss(n,1))));

  }

}

