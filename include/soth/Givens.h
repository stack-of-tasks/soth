#ifndef __SOTH_GIVENS__
#define __SOTH_GIVENS__

#include <Eigen/Core>
#include <Eigen/Jacobi>


namespace soth
{
  using namespace Eigen;

  class Givens
  {
  private:
    typedef PlanarRotation<double> NestedType;

  public:
    Givens();
    Givens(double a, double b, int i, int j, double* z=0);

    void makeGivens(double a, double b, int i, int j, double* z=0);

    // M := M*G.
    template<typename Derived>
    void applyThisOnTheLeft(MatrixBase<Derived> & M) const
    {
      M.applyOnTheRight(i, j, G);
    }

    // M := M*G'.
    template<typename Derived>
    void applyTransposeOnTheLeft(MatrixBase<Derived> & M) const
    {
      M.applyOnTheRight(i, j, G.adjoint());
    }

    // M := G*M.
    template<typename Derived>
    void applyThisOnTheRight(MatrixBase<Derived> & M) const
    {
      M.applyOnTheLeft(i, j, G);
    }

    // M := G'*M.
    template<typename Derived>
    void applyTransposeOnTheRight(MatrixBase<Derived> & M) const
    {
      M.applyOnTheLeft(i, j, G.adjoint());
    }

  private:
    NestedType G;
    int i;
    int j;
  };
}

#endif //__SOTH_GIVENS__