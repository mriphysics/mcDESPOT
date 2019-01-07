/*
 * MEX_mcDESPOT_B0.cpp.
 *
 * Takes input of T1F, T1S, T2F, T2S, kFS, M0F, TRs, Delta, FAs and 
 * outputs 3 FA-by-1 signal vectors. Accounts for B0-effects.
 *
 * This is a MEX file for MATLAB.
*/

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include "signal_SPGR_DE.h"
#include "signal_bSSFP0_B0_DE.h"
#include "signal_bSSFP180_B0_DE.h"
#include "mex.h"

/* The computational routine. */

using namespace Eigen;

void MEX_mcDESPOT(double t1f, double t1s, double t2f, double t2s, double kfs, double m0f, double trspgr, double trssfp, double delta, double *faspgr, double *fassfp, double *SPGR_Sig, double *bSSFP0_Sig, double *bSSFP180_Sig, mwSize n_faspgr, mwSize n_fassfp)

{
  // Assume for modelling.
  int niter = 1;

  SPGR_steady_state_M0 spgr;
  SSFP_steady_state_IV bssfp;
  SSFP_steady_state_180_IV bssfp180;

  spgr.T1_F() = t1f; bssfp.T1_F() = t1f; bssfp180.T1_F() = t1f;
  spgr.T1_S() = t1s; bssfp.T1_S() = t1s; bssfp180.T1_S() = t1s;
  
  bssfp.T2_F() = t2f; bssfp180.T2_F() = t2f;
  bssfp.T2_S() = t2s; bssfp180.T2_S() = t2s;

  spgr.k_FS() = kfs; bssfp.k_FS() = kfs; bssfp180.k_FS() = kfs;  
  spgr.M0_F() = m0f; bssfp.M0_F() = m0f; bssfp180.M0_F() = m0f;

  spgr.TR() = trspgr; bssfp.TR() = trssfp; bssfp180.TR() = trssfp;

  bssfp.Delta() = delta; bssfp180.Delta() = delta;
  
  mwSize i;
  spgr.FA.resize (n_faspgr);
  for (i=0; i<n_faspgr; i++) {
 
    spgr.FA(i) = faspgr[i] * (3.14159265358979323846 / 180.0);
    
  }

  mwSize j;
  bssfp.FA.resize (n_fassfp);
  bssfp180.FA.resize (n_fassfp);
  for (j=0; j<n_fassfp; j++) {
 
    bssfp.FA(j) = fassfp[j] * (3.14159265358979323846 / 180.0);
    bssfp180.FA(j) = fassfp[j] * (3.14159265358979323846 / 180.0);    

  }

  spgr.compute ();
  bssfp.compute();
  bssfp180.compute();
  
  for (int c = 0; c < spgr.FA.size(); c++){
      SPGR_Sig[c] = spgr.Mss_Sig(c);
  }

 for (int c = 0; c < bssfp.FA.size(); c++){
      bSSFP0_Sig[c] = bssfp.Mss_Sig(c);
      bSSFP180_Sig[c] = bssfp180.Mss_Sig(c);
  }

}

/* The gateway function. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double RT1, RT2, RT3, RT4, ER1, FM1, SPGR_TR, SSFP_TR, B0;
    double *SPGR_FA;                                                          
    size_t ncols_spgr;
    double *SSFP_FA;
    size_t ncols_ssfp;                                                           
    double *SPGR, *bSSFP0, *bSSFP180;                            

    /* Check for proper number of arguments. */
    if(nrhs!=11) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","11 inputs required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","3 outputs required.");
    }
    
    /* Get the value of the scalar input.  */
    RT1 = mxGetScalar(prhs[0]);
    RT2 = mxGetScalar(prhs[1]);
    RT3 = mxGetScalar(prhs[2]);
    RT4 = mxGetScalar(prhs[3]);
    ER1 = mxGetScalar(prhs[4]);
    FM1 = mxGetScalar(prhs[5]);
    SPGR_TR = mxGetScalar(prhs[6]);
    SSFP_TR = mxGetScalar(prhs[7]);
    B0 = mxGetScalar(prhs[8]);

    /* Create a pointer to the real data in the input matrices.  */
    SPGR_FA = mxGetPr(prhs[9]);
    SSFP_FA = mxGetPr(prhs[10]);

    /* Get dimensions of the input matrices. */
    ncols_spgr = mxGetN(prhs[9]);
    ncols_ssfp = mxGetN(prhs[10]);

    /* Create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols_spgr,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,(mwSize)ncols_ssfp,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,(mwSize)ncols_ssfp,mxREAL);

    /* Get a pointer to the real data in the output matrix. */
    SPGR = mxGetPr(plhs[0]);
    bSSFP0 = mxGetPr(plhs[1]);
    bSSFP180 = mxGetPr(plhs[2]);

    /* Call the computational routine. */
    MEX_mcDESPOT(RT1,RT2,RT3,RT4,ER1,FM1,SPGR_TR,SSFP_TR,B0,SPGR_FA,SSFP_FA,SPGR,bSSFP0,bSSFP180,(mwSize)ncols_spgr,(mwSize)ncols_ssfp);
    
}
