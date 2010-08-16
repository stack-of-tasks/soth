#include <Eigen/Core>
#include "soth/SubMatrix.hpp"
#include <iostream>
#include "soth/solvers.hpp"

using namespace soth;

void testSolve()
{
  const int n = 6;
  MatrixXd A = MatrixXd::Random(n,n);
  SubMatrix<MatrixXd> P1(A,true,true);
  MatrixXd P2 = A;
  for (int i=0; i<5; ++i)
  {
    int j = rand()%(n);
    int k = rand()%(n);
    P1.permuteCols(j,k);
    P2.col(j).swap(P2.col(k));
  }
  for (int i=0; i<5; ++i)
  {
    int j = rand()%(n);
    int k = rand()%(n);
    P1.permuteRows(j,k);
    P2.row(j).swap(P2.row(k));
  }
  std::cout << (P1 - P2).isZero() << std::endl;
  assert( (P1 - P2).isZero() );

  VectorXd b = VectorXd::Random(n);
  VectorXd b1 = b;
  VectorXd b2 = b;

  soth::solveInPlaceWithLowerTriangular(P1,b1);
  P2.triangularView<Lower>().solveInPlace(b2);
  std::cout << b1.transpose() << std::endl;
  std::cout << b2.transpose() << std::endl;
  std::cout << ((b1-b2).isZero()? "solution is ok": "there is a problem...") << std::endl;
  assert( (b1 - b2).isZero() );
}





int main()
{
  //testSubMatrix();
  //speedTest();

  testSolve();
}
