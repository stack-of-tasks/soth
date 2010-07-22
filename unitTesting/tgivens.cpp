#include "soth/Givens.h"
#include "soth/Algebra.h"

using namespace Eigen;
using namespace soth;

int main()
{
  MatrixXd M = MatrixXd::Random(4,6);
  Givens G1(M(2,0), M(3,0), 2, 3);
  std::cout << (MATLAB)M << std::endl << std::endl;
  G1.applyTransposeOnTheRight(M);
  std::cout << (MATLAB)M << std::endl << std::endl;
  Givens G2(M(1,0), M(2,0), 1, 2);
  G2.applyTransposeOnTheRight(M);
  std::cout << (MATLAB)M << std::endl << std::endl;
  Givens G3(M(0,0), M(1,0), 0, 1);
  G3.applyTransposeOnTheRight(M);
  std::cout << (MATLAB)M << std::endl << std::endl;

  for (int i=M.cols()-1; i>1; --i)
  {
    Givens G(M(0,i-1), M(0,i), i-1, i);
    G.applyThisOnTheLeft(M);
    std::cout << (MATLAB)M << std::endl << std::endl;
  }
}
