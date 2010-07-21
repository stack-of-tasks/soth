#include "soth/BaseY.hpp"

namespace soth
{
  // Empty construction with memory allocation.
  BaseY::BaseY( const unsigned int & insize )
    :isExplicit(false)
    ,size(insize)
    ,matrixExplicit(size,size)
    ,matrixHH()
  {
    matrixHH.reserve(size);
  }

  void BaseY::
  computeExplicitly()
  {
    /* TODO */
    throw "TODO";
  }


  // Y *= Yup. Carefull: there is some recopy here.
  void BaseY::
  composeOnTheRight( const BaseY& Yp )
  {
    /* TODO */
    throw "TODO";
  }
  // Y *= HHup.

  void BaseY::
  composeOnTheRight( const HouseholderSequence & hh )
  {
    if( isExplicit )
      {
	/* TODO */
	throw "TODO";
      }
    else
      {
	matrixHH.append(hh);
      }

  }

}; // namespace soth
