/*
 *  Copyright
 */

#include "soth/Algebra.h"
#include "soth/RotationHouseholder.hpp"
#include <iostream>
#include <Eigen/QR>


int main (int argc, char** argv)
{
  soth::MatrixXd J;
  soth::randMatrix(J,5,9);
  std::cout << "J = " << (soth::MATLAB)J << std::endl;

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Aqr(9,5);
  Aqr.compute(J.transpose());

  soth::MatrixXd Q = Aqr.householderQ();
  std::cout << "Q = " << (soth::MATLAB)Q << std::endl;

  soth::MatrixXd P = Aqr.colsPermutation();
  std::cout << "P = " << (soth::MATLAB)P << std::endl;

  const Eigen::MatrixXd & mQR = Aqr.matrixQR();
  const Eigen::VectorXd & coeff = Aqr.hCoeffs();

  std::cout << "QR = " << (soth::MATLAB)mQR << std::endl;
  std::cout << "b = " << (soth::MATLAB)coeff << std::endl;

  soth::RotationHouseholder hh1;
  hh1.init(mQR,coeff,0 );

  Eigen::VectorXd v = J.row(2);
  hh1.multiplyVector( v );
  hh1.multiplyVector( J.row(2) );
  std::cout << "HJ = " << (soth::MATLAB)v << std::endl;


}
