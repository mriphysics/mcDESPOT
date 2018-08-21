/*
 * MEX_mcDESPOT_SP.cpp.
 *
 * Takes input of T1, T2, M0, TRs, Delta, FAs and 
 * outputs 3 FA-by-1 signal vectors. Accounts for B0-effects.
 *
 * This is a MEX file for MATLAB.
*/

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include "signal_SPGR_SP.h"
#include "signal_bSSFP0_SP.h"
#include "signal_bSSFP180_SP.h"
#include "mex.h"

/* The computational routine. */

using namespace Eigen;

void MEX_mcDESPOT(double t1, double t2, double m0, double trspgr, double trssfp, double delta, double *faspgr, double *fassfp0, double *fassfp180, double *SPGR_Sig, double *bSSFP0_Sig, double *bSSFP180_Sig, mwSize n_faspgr, mwSize n_fassfp0, mwSize n_fassfp180)

{
  // Assume for modelling.
  int niter = 1;

  SPGR_SP_SS spgr;
  SSFP_SP_SS bssfp;
  SSFP_SP_SS_180 bssfp180;

  spgr.T1() = t1; bssfp.T1() = t1; bssfp180.T1() = t1;
  bssfp.T2() = t2; bssfp180.T2() = t2;
  spgr.M0() = m0; bssfp.M0() = m0; bssfp180.M0() = m0;

  spgr.TR() = trspgr; bssfp.TR() = trssfp; bssfp180.TR() = trssfp;
  bssfp.Delta() = delta; bssfp180.Delta() = delta;
  
  mwSize i;
  spgr.FA.resize (n_faspgr);
  for (i=0; i<n_faspgr; i++) {
 
    spgr.FA(i) = faspgr[i] * (3.14159265358979323846 / 180.0);
    
  }

  mwSize j;
  bssfp.FA.resize (n_fassfp0);
  for (j=0; j<n_fassfp0; j++) {
 
    bssfp.FA(j) = fassfp0[j] * (3.14159265358979323846 / 180.0);
    
  }

  mwSize k;
  bssfp180.FA.resize (n_fassfp180);
  for (k=0; k<n_fassfp180; k++) {

    bssfp180.FA(k) = fassfp180[k] * (3.14159265358979323846 / 180.0);

  }  

  spgr.compute ();
  bssfp.compute();
  bssfp180.compute();
  
  for (int c = 0; c < spgr.FA.size(); c++){
    SPGR_Sig[c] = spgr.Mss_Sig(c);
  }

 for (int c = 0; c < bssfp.FA.size(); c++){
    bSSFP0_Sig[c] = bssfp.Mss_Sig(c);
 }

 for (int c = 0; c < bssfp180.FA.size(); c++){
    bSSFP180_Sig[c] = bssfp180.Mss_Sig(c);
 }

}

/* The gateway function. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double RT1, RT2, FM1, SPGR_TR, SSFP_TR, B0;
    double *SPGR_FA;                                                          
    size_t ncols_spgr;
    double *SSFP0_FA;
    size_t ncols_ssfp0;
    double *SSFP180_FA;
    size_t ncols_ssfp180;                                                           
    double *SPGR, *bSSFP0, *bSSFP180;                            

    /* Check for proper number of arguments. */
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","9 inputs required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","3 outputs required.");
    }
    
    /* Get the value of the scalar input.  */
    RT1 = mxGetScalar(prhs[0]);
    RT2 = mxGetScalar(prhs[1]);
    FM1 = mxGetScalar(prhs[2]);
    SPGR_TR = mxGetScalar(prhs[3]);
    SSFP_TR = mxGetScalar(prhs[4]);
    B0 = mxGetScalar(prhs[5]);

    /* Create a pointer to the real data in the input matrices.  */
    SPGR_FA = mxGetPr(prhs[6]);
    SSFP0_FA = mxGetPr(prhs[7]);
    SSFP180_FA = mxGetPr(prhs[8]);

    /* Get dimensions of the input matrices. */
    ncols_spgr = mxGetN(prhs[6]);
    ncols_ssfp0 = mxGetN(prhs[7]);
    ncols_ssfp180 = mxGetN(prhs[8]);

    /* Create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols_spgr,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,(mwSize)ncols_ssfp0,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,(mwSize)ncols_ssfp180,mxREAL);

    /* Get a pointer to the real data in the output matrix. */
    SPGR = mxGetPr(plhs[0]);
    bSSFP0 = mxGetPr(plhs[1]);
    bSSFP180 = mxGetPr(plhs[2]);

    /* Call the computational routine. */
    MEX_mcDESPOT(RT1,RT2,FM1,SPGR_TR,SSFP_TR,B0,SPGR_FA,SSFP0_FA,SSFP180_FA,SPGR,bSSFP0,bSSFP180,(mwSize)ncols_spgr,(mwSize)ncols_ssfp0,(mwSize)ncols_ssfp180);
    
}
