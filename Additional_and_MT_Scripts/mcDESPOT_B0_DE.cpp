#include <iostream>
#include <sstream>
#include <vector>
#include <random>
#include <fstream>
#include "signal_SPGR_DE.h"
#include "signal_bSSFP0_B0_DE.h"
#include "signal_bSSFP180_B0_DE.h"

using namespace Eigen;

inline double as_double (const std::string& value)
{
  double d;
  std::istringstream (value) >> d;
  return d;
}

inline int as_int (const std::string& value)
{
  int d;
  std::istringstream (value) >> d;
  return d;
}


template <typename Derived>
inline void add_noise (DenseBase<Derived>& mat, double sigma)
{
  if (sigma==0.0)
    return;

  static std::random_device rd;
  static std::mt19937 rng (rd());
  static std::normal_distribution<double> normal_dist (0.0, sigma);

  for (int c = 0; c < mat.cols(); c++)
    for (int r = 0; r < mat.rows(); r++)
      mat(r,c) += normal_dist (rng);
}


int main (int argc, char** argv)
{
  double sigma = 0.0;
  int niter = 1;
  
  SPGR_steady_state_M0 spgr;
  SSFP_steady_state_IV bssfp;
  SSFP_steady_state_180_IV bssfp180;
  try {
    for (int n = 1; n < argc; n++) {
      std::string arg (argv[n]);
      auto pos = arg.find ('=');
      if (pos == std::string::npos)
        throw std::string ("invalid argument");
      auto key = arg.substr (0, pos);
      std::transform (key.begin(), key.end(), key.begin(), ::tolower);
      auto value = arg.substr (pos+1);

      if (key == "t1f") {spgr.T1_F() = as_double (value);
                         bssfp.T1_F() = as_double (value);
                         bssfp180.T1_F() = as_double (value);}
      else if (key == "t1s") {spgr.T1_S() = as_double (value);
                              bssfp.T1_S() = as_double (value);
                              bssfp180.T1_S() = as_double (value);}
      else if (key == "m0f") {spgr.M0_F() = as_double (value);
      			      bssfp.M0_F() = as_double (value);
      			      bssfp180.M0_F() = as_double (value);}
      else if (key == "kfs") {spgr.k_FS() = as_double (value);
      			      bssfp.k_FS() = as_double (value);
      			      bssfp180.k_FS() = as_double (value);}
      else if (key == "trspgr") {spgr.TR() = as_double (value);}
      else if (key == "trssfp") {bssfp.TR() = as_double (value);
                                 bssfp180.TR() = as_double (value);}
      else if (key == "t2f") {bssfp.T2_F() = as_double (value);
                              bssfp180.T2_F() = as_double (value);}
      else if (key == "t2s") {bssfp.T2_S() = as_double (value);
                              bssfp180.T2_S() = as_double (value);}
      else if (key == "sigma") sigma = as_double (value);
      else if (key == "n") niter = as_int (value);
      else if (key == "delta") {bssfp.Delta() = as_double (value);
      			        bssfp180.Delta() = as_double (value);}

      else if (key == "faspgr") {
        std::istringstream faspgr_stream (value);
        std::vector<double> faspgr;
        do {
          double d;
          faspgr_stream >> d;
          if (!faspgr_stream.good())
            break;
          faspgr.push_back (d * 3.14159265358979323846 / 180.0);
        } while (faspgr_stream.good());
        spgr.FA.resize (faspgr.size());
        for (int n = 0; n < faspgr.size(); n++)
          spgr.FA[n] = faspgr[n];
      }
      else if (key == "fassfp") {
        std::istringstream fassfp_stream (value);
        std::vector<double> fassfp;
        do {
          double d;
          fassfp_stream >> d;
          if (!fassfp_stream.good())
            break;
          fassfp.push_back (d * 3.14159265358979323846 / 180.0);
        } while (fassfp_stream.good());
        bssfp.FA.resize (fassfp.size());
        bssfp180.FA.resize (fassfp.size());
        for (int n = 0; n < fassfp.size(); n++){
          bssfp.FA[n] = fassfp[n];
          bssfp180.FA[n] = fassfp[n];}
      }
      else
        throw std::string ("unknown key \"" + key + "\" in command-line arguments");

    }
  }
  catch (std::string message) {
    std::cerr << "ERROR: " << message << std::endl;
    return 1;
  }

  MatrixXd sigs_spgr = MatrixXd (spgr.FA.size(), niter);
  for (int n = 0; n < niter; n++) {
    //add_noise (spgr.params, 0.01);
    for (int c = 0; c < spgr.params.size(); c++)
      //if (spgr.params[c] < 0.01)
          //spgr.params[c] = 0.01;
    spgr.compute ();
    sigs_spgr.col(n) = spgr.Mss_Sig;
  }

  MatrixXd sigs_bssfp = MatrixXd (bssfp.FA.size(), niter);
  for (int n = 0; n < niter; n++) {
    //add_noise (bssfp.params, 0.01);
    for (int c = 0; c < bssfp.params.size(); c++)
      //if (bssfp.params[c] < 0.01)
          //bssfp.params[c] = 0.01;
    bssfp.compute ();
    sigs_bssfp.col(n) = bssfp.Mss_Sig;
  }

  MatrixXd sigs_bssfp180 = MatrixXd (bssfp180.FA.size(), niter);
  for (int n = 0; n < niter; n++) {
    //add_noise (bssfp180.params, 0.01);
    for (int c = 0; c < bssfp180.params.size(); c++)
      //if (bssfp180.params[c] < 0.01)
          //bssfp180.params[c] = 0.01;
    bssfp180.compute ();
    sigs_bssfp180.col(n) = bssfp180.Mss_Sig;
  }

  add_noise (sigs_spgr, sigma);
  add_noise (sigs_bssfp, sigma);
  add_noise (sigs_bssfp180, sigma);
 
  std::cerr << "computed " << sigs_spgr.rows() << " x " << sigs_spgr.cols() << std::endl;
  for (int c = 0; c < sigs_spgr.cols(); c++) {
    for (int r = 0; r < sigs_spgr.rows(); r++)
      std::cout << sigs_spgr(r,c) << " ";
    std::cout << std::endl;
  }

  std::cerr << "computed " << sigs_bssfp.rows() << " x " << sigs_bssfp.cols() << std::endl;
  for (int c = 0; c < sigs_bssfp.cols(); c++) {
    for (int r = 0; r < sigs_bssfp.rows(); r++)
      std::cout << sigs_bssfp(r,c) << " ";
    std::cout << std::endl;
  }

  std::cerr << "computed " << sigs_bssfp180.rows() << " x " << sigs_bssfp180.cols() << std::endl;
  for (int c = 0; c < sigs_bssfp180.cols(); c++) {
    for (int r = 0; r < sigs_bssfp180.rows(); r++)
      std::cout << sigs_bssfp180(r,c) << " ";
    std::cout << std::endl;
  }

  std::ofstream file("mcDESPOT_SignalValues.txt");
  if (file.is_open())
  {
    file << sigs_spgr << " ";
    file << std::endl;
    file << sigs_bssfp << " ";
    file << std::endl;
    file << sigs_bssfp180 << " ";
    file << std::endl;
  }
  file.close();

  return 0;
}
