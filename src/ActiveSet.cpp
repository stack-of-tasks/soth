#include "soth/ActiveSet.hpp"

namespace soth
{

  /* ------------------------------------------------------------------------ */
  /* --- CONSTRUCTION ------------------------------------------------------- */
  /* ------------------------------------------------------------------------ */

  ActiveSet::
  ActiveSet( unsigned int nr )
    :v(nr),freerow(nr),freezed(nr),nba(0)
  {
    reset();
  }

  void ActiveSet::
  reset( void )
  {
    std::fill( v.begin(),v.end(),cstref_t(Bound::BOUND_NONE,-1) );
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
    assert( ref<v.size() );
    assert( (v[ref].first == Bound::BOUND_NONE)&&"Constraint has not been properly unactivated." );
    assert( row<v.size() );
    assert( freerow[row] );
    assert( type!=Bound::BOUND_NONE );

    v[ref]=cstref_t(type,row);
    freerow[row]=false; nba++;
  }

  /* Active the given constraint at any free row of ML, and return the
   * position of the selected free row. */
  int ActiveSet::
  activeRow( unsigned int ref, Bound::bound_t type )
  {
    int row = -1;
    /* TODO: it is possible to store the first freeline. */
    for( unsigned int i=0;i<freerow.size();++i )
      {
	if( freerow[i] )
	  { row=i; break; }
      }
    assert( row>=0 );
    active( ref,type,row );
    return row;
  }


  /* Unactive the constraint at line <row> of ML and frees the corresponding line. */
  void ActiveSet::
  unactiveRow( unsigned int row )
  {
    assert( row<v.size() );
    assert( ! freerow[row] );
    for( unsigned int i=0;i<v.size();++i )
      {
	if( v[i].second == row )
	  {
	    assert( v[i].first != Bound::BOUND_NONE );
	    v[i] = cstref_t(Bound::BOUND_NONE,-1);
	  }
      }
    freerow[row]=true;
    nba--;
  }

  /* ------------------------------------------------------------------------ */
  /* --- MISC MODIFIORS ----------------------------------------------------- */
  /* ------------------------------------------------------------------------ */
  /* Pass the constraint to a twin mode. */
  void ActiveSet::
  freeze( unsigned int ref )
  {
    assert( ref<v.size() );
    assert( v[ref].first != Bound::BOUND_NONE );
    assert( v[ref].first != Bound::BOUND_TWIN );
    //v[ref].first = Bound::BOUND_TWIN;
    assert(! freezed[ref] );
    freezed[ref]=true;
  }


  /* --- DEPRECATED --- */
  /*     DPC */void ActiveSet::
  /*     DPC */permuteRows( const VectorXi & P )
    /*   DPC */{
    /*   DPC */  assert(P.size()==nba);
    /*   DPC */  VectorXi Pt(P.size());
    /*   DPC */  for( unsigned int i=0;i<nba;++i ) Pt(P(i))=i;
    /*   DPC */  for( unsigned int i=0;i<v.size();++i )
      /* DPC */    if( isActive(i) ) v[i].second = Pt(v[i].second);
    /*   DPC */}

  /* ------------------------------------------------------------------------ */
  /* --- ACCESSORS ---------------------------------------------------------- */
  /* ------------------------------------------------------------------------ */
  /* Return the number of active constraint. */
  unsigned int ActiveSet::
  nbActive( void ) const { return nba; }

  /* Give the reference of the constraint (ie line in J) located at row <row> of
   * the working space ML. */
  unsigned int ActiveSet::
  whichConstraint( unsigned int r ) const
  {
    for( unsigned int i=0;i<v.size();++i )
      if( v[i].second == r ) return i;

    const bool THE_REQUESTED_ROW_IS_NOT_ACTIVE = false;
    assert( THE_REQUESTED_ROW_IS_NOT_ACTIVE );
  }
  unsigned int ActiveSet::
  where( unsigned int ref ) const
  {
    assert( v[ref].first != Bound::BOUND_NONE );
    return v[ref].second;
  }
  Bound::bound_t ActiveSet::
  whichBound( unsigned int ref ) const
  {
    assert( v[ref].first != Bound::BOUND_NONE );
    return v[ref].first;
  }
  bool ActiveSet::
  isActive( unsigned int ref ) const
  {
    return( v[ref].first != Bound::BOUND_NONE );
  }
  double ActiveSet::
  sign( unsigned int ref ) const
  {
    assert( v[ref].first != Bound::BOUND_NONE );
    assert( v[ref].first != Bound::BOUND_DOUBLE );
    return (v[ref].first==Bound::BOUND_INF)?-1:+1;
  }
  bool ActiveSet::
  isFreezed( unsigned int ref ) const
  {
    assert( ref<v.size() );
    assert( v[ref].first != Bound::BOUND_NONE );
    return  ( v[ref].first == Bound::BOUND_TWIN )||(freezed[ref]);
  }

  /* ------------------------------------------------------------------------ */
  /* --- DISPLAY ------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */


  void ActiveSet::
  disp( std::ostream& os ) const
  {
    for( unsigned int i=0;i<v.size();++ i )
      {
	os << i << ": ";
	if(isActive(i)) os<< v[i].second; else os <<"Unactive";
	os << std::endl;
      }
    for( unsigned int i=0;i<freerow.size();++ i )
      { os << (freerow[i])?"0":"1"; }
    os<<std::endl;
  }


  /* DEPRECATED: Return a compact of the active line, ordered by row values. */
  /*     DPC   */ActiveSet::
  /*     DPC   */operator VectorXi (void) const
    /*   DPC   */{
    /*   DPC   */  if (nba==0)
      /* DPC   */    return VectorXi();
    /*   DPC   */
    /*   DPC   */  VectorXi res(nba);
    /*   DPC   */  int row = 0;
    /*   DPC   */  for( unsigned int i=0;i<v.size();++i )
      /* DPC   */    if(! freerow[i] ) res(row++) = whichConstraint(i);
    /*   DPC   */   // for( unsigned int i=0;i<v.size();++i )
    /*   DPC   */   // 	if( isActive(i) ) res( where(i) ) = i;
    /*   DPC   */  return res;
    /*   DPC   */}


} // namespace soth

