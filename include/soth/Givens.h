#ifndef __SOTH_GIVENS__
#define __SOTH_GIVENS__

#include <Eigen/Core>
#include <Eigen/Jacobi>
#include <vector>


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

    template<typename VectorBase>
    Givens(const VectorBase & v, int i, int j, double* z=0);

  public:
    void makeGivens(double a, double b, int i, int j, double* z=0);

    template<typename VectorBase>
    void makeGivens(const VectorBase & v, int i, int j, double* z=0);

    template<typename VectorBase>
    void makeGivensAndApply(VectorBase & v, int i, int j);

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

    bool applicable( unsigned int size ) const
    {
      bool ci = (i<size);
      bool cj = (j<size);
      assert( (ci&&cj)||((!ci)&&(!cj)) );
      return ci&&cj;
    }

  private:
    NestedType G;
    int i;
    int j;
  };









  class GivensSequence
  {
  public:
    GivensSequence& push(const Givens& g)
    {
      G.push_back(g);
      return *this;
    }


    // M := M*G.
    template<typename Derived>
    void applyThisOnTheLeft(MatrixBase<Derived> & M) const
    {
      for (size_t i=0; i<G.size(); ++i)
        G[i].applyThisOnTheLeft(M);
    }

    // M := M*G.
    template<typename Derived>
    void applyThisOnTheLeftReduced(MatrixBase<Derived> & M) const
    {
      for (size_t i=0; i<G.size(); ++i)
        if( G[i].applicable( M.cols()) ) applyThisOnTheLeft(M);
    }

    // M := M*G'.
    template<typename Derived>
    void applyTransposeOnTheLeft(MatrixBase<Derived> & M) const
    {
      for (size_t i=0; i<G.size(); ++i)
        G[i].applyTransposeOnTheLeft(M);
    }

    // M := G*M.
    template<typename Derived>
    void applyThisOnTheRight(MatrixBase<Derived> & M) const
    {
      for (size_t i=0; i<G.size(); ++i)
        G[i].applyThisOnTheRight(M);
    }

    // M := G'*M.
    template<typename Derived>
    void applyTransposeOnTheRight(MatrixBase<Derived> & M) const
    {
      for (size_t i=0; i<G.size(); ++i)
        G[i].applyTransposeOnTheRight(M);
    }

  private:
    std::vector<Givens> G;
  };


  template<typename VectorBase>
  inline Givens::Givens(const VectorBase & v, int i, int j, double* z)
    :i(i), j(j)
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBase)
    makeGivens(v, i, j, z);
  }

  template<typename VectorBase>
  inline void Givens::makeGivens(const VectorBase & v, int i, int j, double* z)
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBase)
    makeGivens(v(i), v(j), i, j, z);
  }

  template<typename VectorBase>
  inline void Givens::makeGivensAndApply(VectorBase & v, int i, int j)
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBase)
    makeGivens(v(i), v(j), i, j, &v(i));
    v(j) = 0;
  }
}

#endif //__SOTH_GIVENS__
