#include <unsupported/Eigen/MatrixFunctions>
#include "signal_bSSFP180_B0.h"
#include <complex>
#include <iostream>

using namespace Eigen;

void SSFP_steady_state_180_IV::compute ()
{
  // Calculate remaining parameters.
  double R1_S = 1.0/T1_S();
  double R1_F = 1.0/T1_F();
  double R2_S = 1.0/T2_S();
  double R2_F = 1.0/T2_F();

  // For dynamic equilibrium:
  // M0_S = 1 - M0_F; k_SF = (M0_F*k_FS)/M0_S;

  // Define Bloch-McConnell terms:
  VectorXd C(6);
  C << 0.0, 0.0, (R1_F * M0_F()), 0.0, 0.0, (R1_S * M0_S());

  MatrixXd A(6,6);
  A << -R2_F-k_FS(), 0.0, 0.0, k_SF(), 0.0, 0.0,
    0.0, -R2_F-k_FS(), 0.0, 0.0, k_SF(), 0.0,
    0.0, 0.0, -R1_F-k_FS(), 0.0, 0.0, k_SF(),
    k_FS(), 0.0, 0.0, -R2_S-k_SF(), 0.0, 0.0,
    0.0, k_FS(), 0.0, 0.0, -R2_S-k_SF(), 0.0,
    0.0, 0.0, k_FS(), 0.0, 0.0, -R1_S-k_SF();
    
  Mss.resize (FA.size(), 6);
  Mss_Sig.resize (FA.size());

  PartialPivLU<MatrixXd> lu (A);
  VectorXd AinvC(6);
  AinvC = lu.solve (C);

  A *= TR();
  MatrixXd em(6,6);
  em = A.exp();

  const   std::complex<double> i(0.0,1.0); 

  double PCTheta = 3.14159265358979323846 + Delta();

  MatrixXd PCM(6,6);

  PCM << std::cos(PCTheta), -std::sin(PCTheta), 0.0, 0.0, 0.0, 0.0,
  std::sin(PCTheta), std::cos(PCTheta), 0.0, 0.0, 0.0, 0.0, 
  0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, std::cos(PCTheta), -std::sin(PCTheta), 0.0,
  0.0, 0.0, 0.0, std::sin(PCTheta), std::cos(PCTheta), 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

  MatrixXd T(6,6);

  for (int n = 0; n < FA.size(); n++) {

    T << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, std::cos(FA(n)), std::sin(FA(n)), 0.0, 0.0, 0.0,
      0.0, -std::sin(FA(n)), std::cos(FA(n)), 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, std::cos(FA(n)), std::sin(FA(n)),
      0.0, 0.0, 0.0, 0.0, -std::sin(FA(n)), std::cos(FA(n));

    lu.compute (PCM - (em * T));
    Mss.row(n) = lu.solve ((em - MatrixXd::Identity(6,6)) * AinvC);

    // Extract signal component.
    Mss_Sig(n) = std::abs((Mss(n,0) + (i * Mss(n,1))) + (Mss(n,3) + (i * Mss(n,4))));

  }

}

