/*
 *  Copyright
 */

#include "soth/Algebra.h"
#include <iostream>

bool soth::MATLAB::fullPrec = true;

int main (int argc, char** argv)
{
  soth::MatrixXd A(4,5);
  soth::randMatrix(A,9,5);
  std::cout << (soth::MATLAB)A << std::endl;

  //bnu::matrix<double,bnu::column_major> Aqr(9,5);

  soth::bnu::vector<int> E(9);
  soth::bnu::vector<double> b(9);
  boost::numeric::bindings::lapack::geqp(A,E,b);


  soth::bnu::matrix<double,soth::bnu::column_major> Aqr;
  soth::randMatrix(Aqr,9,5);
  std::cout << "Aqr =" << (soth::MATLAB)Aqr << std::endl;

  soth::bnu::vector<double> betas(5); betas.clear();
  soth::bnu::vector<int> orderSe(5); orderSe.clear();
  boost::numeric::bindings::lapack::geqp(Aqr,orderSe,betas);
  std::cout << "qr(Aqr) =" << (soth::MATLAB)Aqr << std::endl;
  std::cout << "b =" << (soth::MATLAB)betas << std::endl;
  std::cout << "E =" << (soth::MATLAB)orderSe << std::endl;
  

}
