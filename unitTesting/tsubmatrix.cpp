#include <Eigen/Core>
namespace Eigen
{
  #include "soth/SubMatrix.h"
}
#include <iostream>
#include <time.h>

using namespace Eigen;

void testSubMatrix()
{
  MatrixXi m(4,5);
  Map<MatrixXi> map(m.data(), m.size(), 1);
  map = VectorXi::LinSpaced(0, m.size()-1, m.size());
  std::cout << "original matrix:" << std::endl << std::endl;
  std::cout << m << std::endl;

  SubMatrix<MatrixXi> p(m);
  std::cout << "creating a permuted matrix with default permutation" << std::endl;
  std::cout << p << std::endl << std::endl;

  std::cout << "permuting row 1 and 2" << std::endl;
  p.permuteRow(1,2);
  std::cout << p << std::endl << std::endl;

  std::cout << "permuting col 2 and 4 then 1 and 4" << std::endl;
  p.permuteCol(2,4);
  p.permuteCol(1,4);
  std::cout << "col permutation indices are " << p.getColIndices().transpose() << " and the permuted matrix is" << std::endl;
  std::cout << p << std::endl << std::endl;

  std::cout << "*****Read access******" << std::endl;
  std::cout << "coef (1,4) is " << p(1,4) << std::endl;
  std::cout << "top right corner of size (2,2) is " << std::endl;
  std::cout << p.topRightCorner(2,2) << std::endl;
  std::cout << "and a triangular lower view" << std::endl;
  MatrixXi pt = p.triangularView<Lower>();
  std::cout << pt << std::endl;


  std::cout << std::endl << std::endl;
  std::cout << "----Now a row permuted only----" << std::endl;
  SubMatrix<MatrixXi, RowPermutation> pr(m);
  std::cout << "creating a row permuted matrix with default permutation" << std::endl;
  std::cout << pr << std::endl << std::endl;

  std::cout << "permuting row 1 and 2" << std::endl;
  pr.permuteRow(1,2);
  std::cout << pr << std::endl << std::endl;

  std::cout << "*****Read access******" << std::endl;
  std::cout << "coef (1,4) is " << pr(1,4) << std::endl;
  std::cout << "top right corner of size (2,2) is " << std::endl;
  std::cout << pr.topRightCorner(2,2) << std::endl;
  std::cout << "and a triangular lower view" << std::endl;
  pt = pr.triangularView<Lower>();
  std::cout << pt << std::endl;

  std::cout << std::endl << std::endl;
  std::cout << "----Now a col submatrix----" << std::endl;
  SubMatrix<MatrixXi, ColPermutation>::ColIndices ind(3);
  ind << 2,4,1;
  SubMatrix<MatrixXi, ColPermutation> pc(m, ind);
  std::cout << "original matrix:" << std::endl << std::endl;
  std::cout << m << std::endl;
  std::cout << "sub matrix" << std::endl;
  std::cout << pc << std::endl << std::endl;
}

void speedTest()
{
  int N=1000;
  int n[] = {5,10,50,100,250};
  for (int i=0; i<5; ++i)
  {
    std::cout << "size " << n[i] << std::endl; 
    MatrixXd A = MatrixXd::Random(n[i],n[i]);
    MatrixXd B = MatrixXd::Random(n[i],n[i]);
    MatrixXd C1(n[i],n[i]);
    MatrixXd C2(n[i],n[i]);
    MatrixXd C3(n[i],n[i]);
    MatrixXd C4(n[i],n[i]);
    SubMatrix<MatrixXd> Prc(A);
    SubMatrix<MatrixXd, RowPermutation> Pr(A);
    SubMatrix<MatrixXd, ColPermutation> Pc(A);

    double dummy;
    clock_t start, stop;
    double total;

    dummy = 0;
    start = clock();
    for (int j=0; j<N; ++j)
    {
      C1.noalias() = A*B;
      dummy += C1(0,0);
    }
    stop = clock();
    total = static_cast<double>(stop-start)/(CLOCKS_PER_SEC*N)*1000;
    std::cout << "normal mult : " << total  << "ms                     dummy=" << dummy << std::endl;

    dummy = 0;
    start = clock();
    for (int j=0; j<N; ++j)
    {
      C2.noalias() = Prc*B;
      dummy += C2(0,0);
    }
    stop = clock();
    total = static_cast<double>(stop-start)/(CLOCKS_PER_SEC*N)*1000;
    std::cout << "Perm*Normal: " << total  << "ms                     dummy=" << dummy << std::endl;

    dummy = 0;
    start = clock();
    for (int j=0; j<N; ++j)
    {
      C3.noalias() = Pr*B;
      dummy += C3(0,0);
    }
    stop = clock();
    total = static_cast<double>(stop-start)/(CLOCKS_PER_SEC*N)*1000;
    std::cout << "RowPerm*Normal : " << total  << "ms                     dummy=" << dummy << std::endl;

    dummy = 0;
    start = clock();
    for (int j=0; j<N; ++j)
    {
      C4.noalias() = Pc*B;
      dummy += C4(0,0);
    }
    stop = clock();
    total = static_cast<double>(stop-start)/(CLOCKS_PER_SEC*N)*1000;
    std::cout << "ColPerm*Normal : " << total  << "ms                     dummy=" << dummy << std::endl;
  }
}

#include "soth/solvers.h"

void testSolve()
{
  const int n = 6;
  MatrixXd A = MatrixXd::Random(n,n);
  SubMatrix<MatrixXd> P1(A);
  MatrixXd P2 = A;
  for (int i=0; i<5; ++i)
  {
    int j = rand()%n;
    int k = rand()%n;
    P1.permuteCol(j,k);
    P2.col(j).swap(P2.col(k));
  }
  for (int i=0; i<5; ++i)
  {
    int j = rand()%n;
    int k = rand()%n;
    P1.permuteRow(j,k);
    P2.row(j).swap(P2.row(k));
  }
  //std::cout << (P1 - P2).isZero() << std::endl;

  VectorXd b = VectorXd::Random(n);
  VectorXd b1 = b;
  VectorXd b2 = b;

  soth::solveInPlaceWithLowerTriangular(P1,b1);
  P2.triangularView<Lower>().solveInPlace(b2);
  std::cout << b1.transpose() << std::endl;
  std::cout << b2.transpose() << std::endl;
  std::cout << ((b1-b2).isZero()? "solution is ok": "there is a problem...") << std::endl;
}



int main()
{
  //testSubMatrix();
  //speedTest();

  testSolve();
}
