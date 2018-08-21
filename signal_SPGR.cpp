#include <unsupported/Eigen/MatrixFunctions>
#include "signal_SPGR.h"
#include <iostream>

using namespace Eigen;

void SPGR_steady_state_M0::compute ()
{

  // Calculate remaining parameters.
  double R1_S = 1.0/T1_S();
  double R1_F = 1.0/T1_F();

  // For dynamic equilibrium:
  // M0_S = 1 - M0_F; k_SF = (M0_F*k_FS)/M0_S;

  // Define Bloch-McConnell terms:
  Vector2d C ((R1_F * M0_F()), (R1_S * M0_S()));

  Matrix2d A;
  A << (-R1_F-k_FS()), k_SF(),
        k_FS(), (-R1_S-k_SF()) ;

  Mss.resize (FA.size(), 2);
  Mss_Sig.resize (FA.size());

  PartialPivLU<Matrix2d> lu (A);
  Vector2d AinvC = lu.solve (C);

  A *= TR();
  Matrix2d em = A.exp();

  for (int n = 0; n < FA.size(); n++) {

    DiagonalMatrix<double,2> T (std::cos(FA(n)), std::cos(FA(n)));

    lu.compute (Matrix2d::Identity() - (em * T));
    Mss.row(n) = lu.solve ((em - Matrix2d::Identity()) * AinvC);

    // Extract signal component.
    Mss_Sig(n) = std::sin(FA(n)) * (Mss(n,0) + Mss(n,1));
  }

}
