#ifndef __SOTH_BASEY__
#define __SOTH_BASEY__


#include "soth/RotationHouseholder.hpp"

namespace soth
{

  class BaseY
  {
  protected:
    typedef MatrixXd::Index Index;

  protected:
    bool isExplicit;
    Index size;
    MatrixXd matrixExplicit;
    MatrixXd householderEssential;
    HouseholderSequence matrixHH;

  public:
    // Empty construction with memory allocation.
    BaseY( const unsigned int & size );

    void computeExplicitly();

  public:
    /* --- Accessor --- */
    MatrixXd& getHouseholderEssential() {return householderEssential;}
    const MatrixXd& getHouseholderEssential() const {return householderEssential;}


    /* --- Multiplier --- */
    /* TODO: when not explicit, it is cheaper to work directly on the
     * result memory, while it is the opposit when explicit. What
     * should we do?
     */

    // v := Y*v = v*Y'.
    template< typename VectorGen >
    void applyThisOnVector( VectorGen & v ) const;
    // v := Y'*v = v*Y.
    template< typename VectorGen >
    void applyTransposeOnVector( VectorGen & v ) const;

    // M := M*Y.
    template< typename Derived >
    void applyThisOnTheLeft( MatrixBase<Derived> & M ) const;
    // M := M*Y'.
    template< typename Derived >
    void applyTransposeOnTheLeft( MatrixBase<Derived> & M ) const;
    // M := Y*M.
    template< typename Derived >
    void applyThisOnTheRight( MatrixBase<Derived> & M ) const;
    // M := Y'*M.
    template< typename Derived >
    void applyTransposeOnTheRight( MatrixBase<Derived> & M ) const;

    // Y *= Yup. Carefull: there is some recopy here.
    void composeOnTheRight( const BaseY& Yp );
    // Y *= HHup.
    void composeOnTheRight( const HouseholderSequence & hh );

  };


  /* --- HEAVY CODE --------------------------------------------------------- */
  /* --- HEAVY CODE --------------------------------------------------------- */
  /* --- HEAVY CODE --------------------------------------------------------- */

  // v := Y*v = v*Y'.
  template< typename VectorGen >
  void BaseY::
  applyThisOnVector( VectorGen & v ) const
  {
    if( isExplicit )
      {
	/*TODO*/ throw "TODO";
      }
    else
      {
	matrixHH.applyThisOnVector(v);
      }
  }

  // v := Y'*v = v*Y.
  template< typename VectorGen >
  void BaseY::
  applyTransposeOnVector( VectorGen & v ) const
  {
    if( isExplicit )
      {
	/*TODO*/ throw "TODO";
      }
    else
      {
	matrixHH.applyTransposeOnVector(v);
      }
  }
  // M := M*Y.
  template< typename Derived >
  void BaseY::
  applyThisOnTheLeft( MatrixBase<Derived> & M ) const
  {
    if( isExplicit )
      {
	/*TODO*/ throw "TODO";
      }
    else
      {
	matrixHH.applyThisOnTheLeft(M);
      }
  }

  // M := M*Y'.
  template< typename Derived >
  void BaseY::
  applyTransposeOnTheLeft( MatrixBase<Derived> & M ) const
  {
    if( isExplicit )
      {
	/*TODO*/ throw "TODO";
      }
    else
      {
	matrixHH.applyTransposeOnTheLeft(M);
      }
  }

  // M := Y*M.
  template< typename Derived >
  void BaseY::
  applyThisOnTheRight( MatrixBase<Derived> & M ) const
  {
    if( isExplicit )
      {
	/*TODO*/ throw "TODO";
      }
    else
      {
	matrixHH.applyThisOnTheRight(M);
      }
  }

  // M := Y'*M.
  template< typename Derived >
  void BaseY::
  applyTransposeOnTheRight( MatrixBase<Derived> & M ) const
  {
    if( isExplicit )
      {
	/*TODO*/ throw "TODO";
      }
    else
      {
	matrixHH.applyTransposeOnTheRight(M);
      }
  }


};




#endif //#ifndef __SOTH_BASEY__
