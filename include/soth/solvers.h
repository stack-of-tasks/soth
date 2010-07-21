#include <Eigen/Core>

using namespace Eigen;

namespace soth
{
  template<typename Derived1, typename Derived2>
  inline void solveInPlaceWithLowerTriangular(const MatrixBase<Derived1>& lhs, MatrixBase<Derived2>& rhs)
  {
    //forward substitution, colum version
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(MatrixBase<Derived2>)
    assert(lhs.rows() == lhs.cols());
    assert(rhs.size() == lhs.rows());
    const int n = lhs.rows();
    for (int i=0; i<n-1; ++i)
    {
      rhs[i] /= lhs(i,i);
      rhs.tail(n-i-1) = rhs.tail(n-i-1) - rhs[i]* lhs.col(i).tail(n-i-1);
    }
    rhs[n-1] /= lhs(n-1,n-1);
  }
}