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
  std::cout << (soth::MATLAB)A << std::endl;

  //Eigen::QR<Eigen::MatrixXd> Aqr(A);
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Aqr(A);

}
