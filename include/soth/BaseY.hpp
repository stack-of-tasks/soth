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
    MatrixXd matrixExplicit;
    HouseholderSequence matrixHH;
    Index size;

  public:
    // Empty construction with memory allocation.
    BaseY( const unsigned int & size );

    void computeExplicitly();

    // a := Y*a
    template< typename VectorGen >
    void applyThisOnTheRight( VectorGen & a );
    // a := Y'*a = a'*Y
    template< typename VectorGen >
    void applyThisOnTheLeft( VectorGen & a );

    // Y *= Yup. Carefull: there is some recopy here.
    void composeOnTheRight( const BaseY& Yp );
    // Y *= HHup.
    void composeOnTheRight( const HouseholderSequence & hh );

  };

};




#endif //#ifndef __SOTH_BASEY__
