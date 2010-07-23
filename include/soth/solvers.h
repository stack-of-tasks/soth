#include <Eigen/Core>

using namespace Eigen;

namespace soth
{

  /* Solve x in the problem Ax=b, for generic matrix and vector A and b, with A=<lhs> and b=<rhs>, and store the result
   * in v. A is considered to be upper triangular and full rank. 
   * Its strictly upper part is considered to be 0 and the solver performs a forward substitution.
   */
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
      rhs.tail(n-i-1) -= rhs[i]* lhs.col(i).tail(n-i-1);
    }
    rhs[n-1] /= lhs(n-1,n-1);
  }


  /* Solve x in the problem Ax=b, for generic matrix and vector A and b, with A=<lhs> and b=<rhs>, and store the result
   * in v. A is considered to be upper triangular and full rank. 
   * Its strictly lower part is considered to be 0 and the solver performs a backward substitution.
   */
  template<typename Derived1, typename Derived2>
  inline void solveInPlaceWithUpperTriangular(const MatrixBase<Derived1>& lhs, MatrixBase<Derived2>& rhs)
  {
    //backward substitution, colum version
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(MatrixBase<Derived2>)
    assert(lhs.rows() == lhs.cols());
    assert(rhs.size() == lhs.rows());
    const int n = lhs.rows();
    for (int i=n-1; i>0; --i)
    {
      rhs[i] /= lhs(i,i);
      rhs.head(i) -= rhs[i]* lhs.col(i).head(i);
      std::cout << rhs.transpose() << std::endl;
    }
    rhs[0] /= lhs(0,0);
  }
}
