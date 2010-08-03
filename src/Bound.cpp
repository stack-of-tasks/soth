#include <assert.h>
#include <cmath>
#include "soth/Bound.hpp"

namespace soth
{
  Bound::Bound( void )
    : type(BOUND_NONE),valInf(0),valSup(0),valTwin(valInf)
  {
  }

  Bound::Bound( const Bound & clone )
    : type(clone.type),valInf(clone.valInf),valSup(clone.valSup),valTwin(valInf)
  {
  }

  Bound::Bound( const double & val, bound_t inType )
    : type(inType),valInf(val),valSup(val),valTwin(valInf)
  {
    assert( inType != BOUND_DOUBLE );
    assert( inType != BOUND_NONE );

    type = inType;
  }

  Bound::Bound( const double & inValInf, const double & inValSup )
    : type(BOUND_DOUBLE),valInf(inValInf),valSup(inValSup),valTwin(valInf)
  {
  }

  const double& Bound::
  getBound( bound_t inType ) const
  {
    assert( inType != BOUND_NONE );
    assert( inType != BOUND_DOUBLE );

    switch( inType )
      {
      case BOUND_INF:
	assert( (type==BOUND_INF)||(type==BOUND_DOUBLE) );
	return valInf;
      case BOUND_SUP:
	assert( (type==BOUND_SUP)||(type==BOUND_DOUBLE) );
	return valSup;
      case BOUND_TWIN:
	assert(type==BOUND_TWIN);
	return valTwin;
      }
  }

  /* Return the bound that is violated, NONE if bound are OK.
   * In case of twin-bounds, no check is performed, NONE is always returned. */
  Bound::bound_t Bound::
  check( const double & val ) const
  {
    assert( type!=BOUND_NONE );
    switch( type )
      {
      case BOUND_INF:
	if( val<valInf ) return BOUND_INF;
	break;
      case BOUND_SUP:
	if( valSup<val ) return BOUND_SUP;
	break;
      case BOUND_TWIN:
	break;
      case BOUND_DOUBLE:
	if( val<valInf ) return BOUND_INF;
	if( valSup<val ) return BOUND_SUP;
      }
    return BOUND_NONE;
  }

  /* Return the bound that is violated, NONE if bound are OK.
   * In case of twin-bounds, no check is performed, NONE is always returned. */
  Bound::bound_t Bound::
  checkSaturation( const double & val,const double & EPSILON ) const
  {
    assert( type!=BOUND_NONE );
    switch( type )
      {
      case BOUND_INF:
	if( std::abs(val-valInf)<EPSILON ) return BOUND_INF;
	break;
      case BOUND_SUP:
	if( std::abs(val-valSup)<EPSILON ) return BOUND_SUP;
	break;
      case BOUND_TWIN:
	break;
      case BOUND_DOUBLE:
	if( std::abs(val-valInf)<EPSILON ) return BOUND_INF;
	if( std::abs(val-valSup)<EPSILON ) return BOUND_SUP;
	break;
      }
    return BOUND_NONE;
  }

  Bound& Bound::operator= ( const Bound& clone )
  {
    valInf=clone.valInf;
    valSup=clone.valSup;
    type=clone.type;
    return *this;
  }

  Bound& Bound::operator= ( const double & val)
  {
    valTwin = val;
    type= BOUND_TWIN;
    return *this;
  }

  Bound& Bound::operator= ( const std::pair<double,double> & val)
  {
    valInf = val.first; valSup = val.second;
    type= BOUND_DOUBLE;
    return *this;
  }


  std::ostream& operator<< (std::ostream& os, const Bound&b )
  {
    switch( b.type )
      {
      case Bound::BOUND_INF:
	os << "[ " << b.valInf << ", +inf ]";
	break;
      case Bound::BOUND_SUP:
	os << "[ -inf ," << b.valSup << " ]";
	break;
      case Bound::BOUND_TWIN:
	os << "[ " << b.valTwin <<", " << b.valTwin << " ]";
	break;
      case Bound::BOUND_DOUBLE:
	os << "[ " << b.valInf <<", " << b.valSup << " ]";
	break;
      case Bound::BOUND_NONE:
	os << "[ -inf, +inf ]";
      }
    return os;
  }
  std::ostream& operator<< (std::ostream& os, const bound_vector_t& t)
  {
    os << "[ ";
    for( bound_vector_t::const_iterator iter=t.begin();
	 iter!=t.end();++iter )
      {
	os << *iter;
	if( iter+1 != t.end() ) os <<"; "; else os<<" ]";
      }
    return os;
  }



}; // namespace soth
