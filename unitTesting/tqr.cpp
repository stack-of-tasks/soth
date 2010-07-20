/*
 *  Copyright
 */

#include "soth/Algebra.h"
#include <iostream>
#include <Eigen/QR>


int main (int argc, char** argv)
{
  soth::MatrixXd A(4,5);
  soth::randMatrix(A,9,5);
  std::cout << "A = " << (soth::MATLAB)A << std::endl;

  //Eigen::QR<Eigen::MatrixXd> Aqr(A);
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Aqr(5,4);

  Aqr.compute(A.transpose());
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Aqr2(A);

  soth::MatrixXd Q = Aqr.householderQ();
  std::cout << "Q = " << (soth::MATLAB)Q << std::endl;

  soth::MatrixXd P = Aqr.colsPermutation();
  std::cout << "P = " << (soth::MATLAB)P << std::endl;

}
