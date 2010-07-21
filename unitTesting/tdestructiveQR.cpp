#include <Eigen/Core>
#include <Eigen/Householder>
#include <Eigen/QR>
namespace Eigen
{
  #include "soth/DestructiveColPivQR.h"
}
#include <iostream>

using namespace Eigen;

void testColPivForCod()
{
  const int n=7;
  const int m=6;
  MatrixXd M = MatrixXd::Random(m,n);
  MatrixXd Y = MatrixXd::Zero(n,n);

  std::cout << "initial matrix M: " << std::endl;
  std::cout << M << std::endl << std::endl << std::endl;

  std::cout << "***** Reference Decomposition of M' *****" << std::endl;
  ColPivHouseholderQR<MatrixXd> ref(M.transpose());
  MatrixXd resR = ref.matrixQR().triangularView<Upper>();
  MatrixXd resQ = ref.matrixQR().triangularView<StrictlyLower>();
  std::cout << "R_ref = " << std::endl;
  std::cout << resR << std::endl << std::endl;
  std::cout << "householder essential H_ref= " << std::endl;
  std::cout << resQ << std::endl << std::endl;

  
  std::cout << "***** Tested Decomposition of M' *****" << std::endl;
  DestructiveColPivQR<Transpose<MatrixXd>, Block<MatrixXd> > qr(M.transpose(),Y.block(0,0,n,m));
  std::cout << "R (in place, no transposition) = " << std::endl;
  std::cout << M.transpose() << std::endl<<std::endl;
  std::cout << "permuted R = " << std::endl;
  std::cout << M.transpose()*qr.colsPermutation() << std::endl<<std::endl;
  std::cout << "householder essential H= " << std::endl;
  std::cout << Y << std::endl<<std::endl;
  std::cout << "QR = " << std::endl;
  std::cout << qr.householderQ()*M.transpose() << std::endl << std::endl;

  std::cout  << std::endl << "***** Check *****" << std::endl;
  std::cout << "correct R (R_ref*P_ref' = R): " 
            << ((resR*ref.colsPermutation().transpose()-M.transpose()).isZero()? "true" : "false") << std::endl;
  std::cout << "correct householder (Href = H): " 
            << ((resQ-Y.block(0,0,n,m)).isZero()? "true" : "false") << std::endl;
  std::cout << "correct permutations (P_ref*P = I): " 
            << ((ref.colsPermutation().transpose()*qr.colsPermutation()).toDenseMatrix().isIdentity()? "true" : "false") << std::endl;
}


void main()
{
  testColPivForCod();
}