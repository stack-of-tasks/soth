#include "soth/ActiveSet.hpp"

namespace soth
{

  /* ------------------------------------------------------------------------ */
  /* --- CONSTRUCTION ------------------------------------------------------- */
  /* ------------------------------------------------------------------------ */

  ActiveSet::
  ActiveSet( unsigned int nr )
    :cstMap(nr),cstMapInv(nr),freerow(nr),freezed(nr),nba(0)
  {
    reset();
  }

  void ActiveSet::
  reset( void )
  {
    std::fill( cstMap.begin(),cstMap.end(),cstref_t(Bound::BOUND_NONE,-1) );
    std::fill( cstMapInv.begin(),cstMapInv.end(),size() );
    std::fill( freerow.begin(),freerow.end(),true );
    std::fill( freezed.begin(),freezed.end(),false );
    nba=0;
  }

  /* ------------------------------------------------------------------------ */
  /* --- ACTIVE ------------------------------------------------------------- */
  /* ------------------------------------------------------------------------ */

  /* Active the constraint <ref> (ie in J(ref,:)) at line <row> (ie in ML(row,:))
   * with bound type. */
  void ActiveSet::
  active( unsigned int ref, Bound::bound_t type, unsigned int row )
  {
    assert( ref<size() );
    assert( (cstMap[ref].first == Bound::BOUND_NONE)&&"Constraint has not been properly unactivated." );
    assert( row<size() );
    assert( freerow[row] );
    assert( cstMapInv[row]==size() );
    assert( type!=Bound::BOUND_NONE );

    cstMap[ref]=cstref_t(type,row);
    cstMapInv[row]=ref;
    freerow[row]=false; nba++;
  }

  /* Active the given constraint at any free row of ML, and return the
   * position of the selected free row. */
  unsigned int ActiveSet::
  activeRow( unsigned int ref, Bound::bound_t type )
  {
    const unsigned int row = getAFreeRow();
    assert( (row>=0)&&(row<size()) );
    active( ref,type,row );
    return row;
  }

  /* Unactive the constraint at line <row> of ML and frees the corresponding line. */
  void ActiveSet::
  unactiveRow( unsigned int row )
  {
    assert( row<size() );

    const unsigned int & cst = cstMapInv[row];
    assert( cst<size() );
    assert( cstMap[cst].first != Bound::BOUND_NONE );

    cstMap[cst] = cstref_t(Bound::BOUND_NONE,-1);
    freeARow(row);
    cstMapInv[row] = size();
    nba--;
  }

  /* Pass the constraint to a twin mode. */
  void ActiveSet::
  freeze( unsigned int ref )
  {
    assert( ref<size() );
    assert( cstMap[ref].first != Bound::BOUND_NONE );
    assert( cstMap[ref].first != Bound::BOUND_TWIN );
    assert(! freezed[ref] );

    freezed[ref]=true;
  }


  bool ActiveSet::
  isRowFree( unsigned int row )
  {
    return freerow[row];
  }
  unsigned int ActiveSet::
  getAFreeRow( void )
  {
    /* TODO: it is possible to store the first freeline. */
    assert( nba < size() );
    for( unsigned int row=0;row<size();++row )
      {
	if( freerow[row] ) return row;
      }
  }
  void ActiveSet::
  freeARow( unsigned int row )
  {
    assert( ! freerow[row] );
    freerow[row]=true;
  }



  /* ------------------------------------------------------------------------ */
  /* --- ACCESSORS ---------------------------------------------------------- */
  /* ------------------------------------------------------------------------ */

  /* Give the reference of the constraint (ie line in J) located at row <row> of
   * the working space ML. */
  unsigned int ActiveSet::
  mapInv( unsigned int row ) const
  {
    assert( row<size() );
    assert( (cstMapInv[row]<size())&&"THE_REQUESTED_ROW_IS_NOT_ACTIVE" );
    return cstMapInv[row];
  }
  unsigned int ActiveSet::
  map( unsigned int ref ) const
  {
    assert( isActive(ref) );
    return cstMap[ref].second;
  }
  Bound::bound_t ActiveSet::
  whichBound( unsigned int ref,bool checkActive ) const
  {
    assert( isActive(ref) );
    const Bound::bound_t & res = cstMap[ref].first;
    if( checkActive ) { assert( (res!=Bound::BOUND_NONE)&&(res!=Bound::BOUND_DOUBLE) ); }
    return res;
  }
  bool ActiveSet::
  isActive( unsigned int ref ) const
  {
    assert( ref<size() );
    return( cstMap[ref].first != Bound::BOUND_NONE );
  }
  double ActiveSet::
  sign( unsigned int ref ) const
  {
    assert( isActive(ref) );
    assert( cstMap[ref].first != Bound::BOUND_DOUBLE );
    return (cstMap[ref].first==Bound::BOUND_INF)?-1:+1;
  }
  bool ActiveSet::
  isFreezed( unsigned int ref ) const
  {
    assert( isActive(ref) );
    return ( cstMap[ref].first == Bound::BOUND_TWIN )||(freezed[ref]);
  }
  /*: Return a compact of the active line, ordered by row values. */
  VectorXi ActiveSet::
  getIndirection(void) const
  {
    if (nba==0)
      return VectorXi();

    VectorXi res(nba);
    int row = 0;
    for( unsigned int i=0;i<size();++i )
      if(! freerow[i] ) res(row++) = whichConstraint(i);
    return res;
  }


  /* ------------------------------------------------------------------------ */
  /* --- DISPLAY ------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */


  void ActiveSet::
  disp( std::ostream& os, bool classic ) const
  {
    if( classic )
      {
	os << " [ ";
	for( unsigned int r=0;r<size();++r )
	  {
	    if( freerow[r] ) continue;
	    const unsigned int cst = mapInv(r);
	    if( whichBound(cst) == Bound::BOUND_INF ) os << "-";
	    else if( whichBound(cst) == Bound::BOUND_SUP ) os << "+";
	    os <<r << " ";
	  }
	os << " ]";
      }
    else
      {
	for( unsigned int i=0;i<size();++ i )
	  {
	    os << i << ": ";
	    if(isActive(i)) os<< cstMap[i].second; else os <<"Unactive";
	    os << std::endl;
	  }
	for( unsigned int i=0;i<freerow.size();++ i )
	  { os << (freerow[i])?"0":"1"; }
	os<<std::endl;
      }
  }
  std::ostream& operator<< ( std::ostream & os,const ActiveSet& as )
  {    as.disp(os); return os; }

  /* ------------------------------------------------------------------------ */
  /* --- DEPRECATED --------------------------------------------------------- */
  /* ------------------------------------------------------------------------ */

  /* --- DEPRECATED --- */
  /*     DPC */void ActiveSet::
  /*     DPC */permuteRows( const VectorXi & P )
    /*   DPC */{
    /*   DPC */  assert( false&&"DEPRECATED" );
    /*   DPC */  assert(P.size()==nba);
    /*   DPC */  VectorXi Pt(P.size());
    /*   DPC */  for( unsigned int i=0;i<nba;++i ) Pt(P(i))=i;
    /*   DPC */  for( unsigned int i=0;i<size();++i )
      /* DPC */    if( isActive(i) ) cstMap[i].second = Pt(cstMap[i].second);
    /*   DPC */}

 
} // namespace soth

