#include "soth/BasicStage.hpp"

namespace soth
{


  std::ostream& operator<< (std::ostream& os, const VectorBound& t)
  {
    os << "[ ";
    typedef VectorBound::Index Index;
    EI_FOREACH( i,t  )
      {
	os << t[i];
	if( i+1<t.size() ) os <<"; "; else os<<" ]";
      }
    return os;
  }

  /* --- STAGE -------------------------------------------------------------- */

  BasicStage::
  BasicStage( const unsigned int innr, const unsigned int innc,
	      const double * Jdata, const Bound * bdata, const BaseY& Y )
    :boundsInternal(0)
    ,Jmap( Jdata,innr,innc )
    ,boundsMap( bdata,innr )

    ,J( Jmap ), bounds( boundsMap )
    ,nr(innr),nc(innc)

    ,Y(Y)
  {}

  BasicStage::
  BasicStage( const unsigned int innr, const unsigned int innc,
	      const double * Jdata, const BaseY& Y )
    :boundsInternal(innr)
    ,Jmap( Jdata,innr,innc )
    ,boundsMap( boundsInternal.data(),innr )

    ,J( Jmap ), bounds( boundsMap )
    ,nr(innr),nc(innc)

    ,Y(Y)
  {}

  BasicStage::
  BasicStage( const MatrixXd & inJ, const VectorBound & inbounds, const BaseY& inY  )
    :boundsInternal(0)
    ,Jmap( inJ.data(),inJ.rows(),inJ.cols() )
    ,boundsMap( inbounds.data(),inbounds.size(),1)

    ,J( Jmap ), bounds( boundsMap )
    ,nr( inJ.rows() ), nc( inJ.cols() )

    ,Y(inY)
  {
    assert( inbounds.size()==int(nr) );
  }

  void BasicStage::
  set( const MatrixXd & inJ, const VectorBound & inbounds )
  {
    assert( inJ.rows() == (int)nr && inJ.cols() == (int)nc );
    assert( inbounds.size() == (int)nr );
    set( inJ.data(),inbounds.data() );
  }

  void BasicStage::
  set( const double * Jdata, const Bound * bdata )
  {
    new (&Jmap) MapXd( Jdata,nr,nc );
    new (&boundsMap) MapBound( bdata,nr );
  }

  MatrixXd BasicStage::
  getJ( void )
  {
    return J;
  }

  VectorBound BasicStage::
  getBounds( void )
  {
    return bounds;
  }

  VectorBound& BasicStage::
  getBoundsInternal()
  {
    assert( boundsInternal.data() == boundsMap.data() );
    return boundsInternal;
  }

} // namespace soth
