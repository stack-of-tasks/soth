#ifndef __SOTH_ASET__
#define __SOTH_ASET__

#include "soth/Bound.hpp"

namespace soth
{

 class ActiveSet
  {
  public:
    ActiveSet( unsigned int nr )
      :v(nr),freerow(nr),nba(0)
    {
      std::fill( v.begin(),v.end(),cstref_t(Bound::BOUND_NONE,-1) );
      std::fill( freerow.begin(),freerow.end(),true );
    }

    /* Return the number of active constraint. */
    unsigned int nbActive( void ) const { return nba; }

    /* Active the constraint <ref> (ie in J(ref,:)) at line <row> (ie in ML(row,:))
     * with bound type. */
    void active( unsigned int ref, Bound::bound_t type, unsigned int row )
    {
      assert( ref<v.size() );
      assert( v[ref].first == Bound::BOUND_NONE );
      assert( row<v.size() );
      assert( freerow[row] );
      assert( type!=Bound::BOUND_NONE );

      v[ref]=cstref_t(type,row);
      freerow[row]=false; nba++;
    }

    /* Active the given constraint at any free row of ML, and return the
     * position of the selected free row. */
    int activeRow( unsigned int ref, Bound::bound_t type )
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
    void unactiveRow( unsigned int row )
    {
      assert( row<v.size() );
      assert( ! freerow[row] );
      for( unsigned int i=0;i<v.size();++i )
	{
	  // if( v[i].second>row ) v[i].second--;
	  // else
	  if( v[i].second == row )
	    {
	      assert( v[i].first != Bound::BOUND_NONE );
	      v[i] = cstref_t(Bound::BOUND_NONE,-1);
	    }
	}
      freerow[row]=true;
      nba--;
    }

    void permuteRows( const VectorXi & P )
    {
      assert(P.size()==nba);
      VectorXi Pt(P.size());
      for( unsigned int i=0;i<nba;++i ) Pt(P(i))=i;
      for( unsigned int i=0;i<v.size();++i )
	if( isActive(i) ) v[i].second = Pt(v[i].second);
    }

    /* Give the reference of the constraint (ie line in J) located at row <row> of
     * the working space ML. */
    unsigned int whichConstraint( unsigned int r ) const
    {
      for( unsigned int i=0;i<v.size();++i )
	if( v[i].second == r ) return i;

      const bool THE_REQUESTED_ROW_IS_NOT_ACTIVE = false;
      assert( THE_REQUESTED_ROW_IS_NOT_ACTIVE );
    }
    unsigned int where( unsigned int ref ) const
    {
      assert( v[ref].first != Bound::BOUND_NONE );
      return v[ref].second;
    }
    Bound::bound_t whichBound( unsigned int ref ) const
    {
      assert( v[ref].first != Bound::BOUND_NONE );
      return v[ref].first;
    }
    bool isActive( unsigned int ref ) const
    {
      return( v[ref].first != Bound::BOUND_NONE );
    }


    /* Return a compact of the active line, ordered by row values. */
    operator VectorXi (void) const
    {
      if (nba==0)
        return VectorXi();

      VectorXi res(nba);
      int row = 0;
      for( unsigned int i=0;i<v.size();++i )
	if(! freerow[i] ) res(row++) = whichConstraint(i);
      // for( unsigned int i=0;i<v.size();++i )
      // 	if( isActive(i) ) res( where(i) ) = i;
      return res;
    }

    void disp( std::ostream& os ) const
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

  protected:
    typedef std::pair<Bound::bound_t,int> cstref_t;
    typedef std::vector<cstref_t> cstref_vector_t;
    cstref_vector_t v;
    std::vector<bool> freerow;
    unsigned int nba;
  };


} // namespace soth

#endif
