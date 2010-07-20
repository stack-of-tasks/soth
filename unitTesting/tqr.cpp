/*
 *  Copyright
 */

#include "soth/Algebra.h"
#include "soth/RotationHouseholder.hpp"
#include <iostream>
#include <Eigen/QR>

template< typename Derived >
void toto( Eigen::MatrixBase<Derived> & v ) {  v(1) = 1; }

template< typename VG>
void tata( VG& v ) {  v(1) = 1; }



int main (int argc, char** argv)
{
  Eigen::MatrixXd J;
  soth::randMatrix(J,5,9);
  std::cout << "J = " << (soth::MATLAB)J << std::endl;

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Aqr(9,5);
  Aqr.compute(J.transpose());

  std::cout << "Q = " << (soth::MATLAB)(soth::MatrixXd) Aqr.householderQ()  << std::endl;
  std::cout << "P = " << (soth::MATLAB)(soth::MatrixXd) Aqr.colsPermutation() << std::endl;

  const Eigen::MatrixXd & mQR = Aqr.matrixQR();
  const Eigen::VectorXd & coeff = Aqr.hCoeffs();

  std::cout << "QR = " << (soth::MATLAB)mQR << std::endl;
  std::cout << "b = " << (soth::MATLAB)coeff << std::endl;


  Eigen::MatrixXd Jt = J.transpose();

  soth::HouseholderSequence hh( mQR,coeff );
  hh.applyThisOnTheLeft(J);
  hh.applyTransposeOnTheRight(Jt);


  std::cout << "JQ = " << (soth::MATLAB)J << std::endl;
  std::cout << "QtJt = " << (soth::MATLAB)Jt << std::endl;


}
