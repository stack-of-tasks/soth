#ifndef __SOTH_BASEY__
#define __SOTH_BASEY__

namespace soth
{
  
  class BaseY
  {
  protected:
    bool isExplicit;
    MatrixXd matrixExplicit;
    householder_sequence_t matrixHH;

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

    // Y *= Yup
    void composeOnTheRight( const BaseY& Yp );


  };





};




#endif //#ifndef __SOTH_BASEY__
