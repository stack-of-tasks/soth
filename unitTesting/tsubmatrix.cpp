#include <Eigen/Core>
namespace Eigen
{
  #include "soth/SubMatrix.h"
}
#include <iostream>

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


int main()
{
  testSubMatrix();
}