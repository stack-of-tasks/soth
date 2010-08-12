#ifndef __SOTH_ASET__
#define __SOTH_ASET__

#include "soth/Algebra.h"
#include "soth/Bound.hpp"

namespace soth
{

  /* ActiveSet is a invertible map that gives the row where an active constraint is stored,
   * ie MAP such as J(MAP,:) == WMLY. The map is invertible, which means that it is possible
   * to acces to the constraint reference of a given row of WMLY. It also stores the type of the
   * active constraint (+/-/=), and the free/occupied lines of WMLY.
   * Take care: WMLY is supposed to be an indirect matrix build on W_MLY, in which case
   * the active set is not aware of the indirection (ie: to access to the constraint of WMLY(row),
   * mapInv should be call with W.indirectRow(row), and not directly with row). Similarly, the
   * row where an active constraint is stored in WMLY will need to inverse the W.indirectRow map, which
   * is not done in the class.
   */
  class ActiveSet
  {
  public: /* --- Construction --- */
    ActiveSet( unsigned int nr );
    void                  reset( void );

  public: /* --- Active set management --- */
    /* Active the constraint <ref> (ie in J(ref,:)) at line <row> (ie in ML(row,:))
     * with bound type. */
    void                  active( unsigned int ref, Bound::bound_t type, unsigned int row );
    /* Active the given constraint at any free row of ML, and return the
     * position of the selected free row. */
    unsigned int          activeRow( unsigned int ref, Bound::bound_t type );
    /* Unactive the constraint at line <row> of ML and frees the corresponding line. */
    void                  unactiveRow( unsigned int row );
    /* Pass the constraint to a twin mode. */
    void                  freeze( unsigned int ref );


  public: /* --- Accessors --- */
    /* Return the number of active constraint. */
    inline unsigned int   nbActive( void ) const { return nba; }
    inline unsigned int   size( void ) const { return cstMap.size(); }
    bool                  isFreezed( unsigned int ref ) const;
    bool                  isActive( unsigned int ref ) const;
    Bound::bound_t        whichBound( unsigned int ref,bool checkActive=false ) const;
    double                sign( unsigned int ref ) const;
    /* Map: from the cst ref, give the row of J.
     * Map inversion: give the reference of the constraint (ie line in J) located at row <row> of
     * the working space ML. */
    unsigned int          map( unsigned int ref ) const;
    unsigned int          mapInv( unsigned int row ) const;
    VectorXi              getIndirection( void ) const;

    /* For compatibility */
    inline unsigned int   where( unsigned int ref ) const { return map(ref); }
    inline unsigned int   whichConstraint( unsigned int row ) const { return mapInv(row); }

  public: /* --- Display --- */
    /* Return a compact of the active line, ordered by row values. */
    void                  disp( std::ostream& os,bool classic=true ) const;

  public: /* --- Deprecated --- */
    /* DEPRECATED*/void   permuteRows( const VectorXi & P );

  protected:
    typedef std::pair<Bound::bound_t,int> cstref_t;
    typedef std::vector<cstref_t> cstref_vector_t;
    typedef std::vector<unsigned int> mapinv_vector_t;
  protected:
    cstref_vector_t cstMap;
    mapinv_vector_t cstMapInv;
    std::vector<bool> freerow,freezed;
    unsigned int nba;

  protected: /* Internal management */
    bool isRowFree( unsigned int row );
    unsigned int getAFreeRow( void );
    void freeARow( unsigned int row );

  public:
    inline operator VectorXi (void) const {  return getIndirection(); }
    friend std::ostream& operator<< ( std::ostream & os,const ActiveSet& as );
  };


  /* The previous class is not aware of the indirection built upon WMLY. This indirection
   * is added in this derivation, to make it transparent to the user.
   */
  template< typename AS,typename Indirect >
  class SubActiveSet
    : public AS
  {
  public:
    SubActiveSet( unsigned int nr );
    SubActiveSet( unsigned int nr, Indirect& idx );
    SubActiveSet( const SubActiveSet& clone );

    inline bool           ownIndirection(void) const { return &self_indirect == &indirect; }

  public:
    void                  reset( void );
    unsigned int          activeRow( unsigned int ref, Bound::bound_t type );
    void                  unactiveRow( unsigned int row );
    unsigned int          mapInv( unsigned int row ) const;
    unsigned int          map( unsigned int ref ) const;
    /* For compatibility */
    inline unsigned int   whichConstraint( unsigned int row ) const { return mapInv(row); }
    void                  disp( std::ostream& os, bool classic=true ) const;

    void                  setInitialActivation( const AS& as0 );

  public:
    using AS::            size;
    using AS::            whichBound;

  protected: /* Protected for clarity. */
    void                  active( unsigned int ref, Bound::bound_t type, unsigned int row );
    unsigned int          pushIndirectBack( unsigned int rowup );

  protected:
    Indirect self_indirect;
    Indirect& indirect;
    bool isEmpty;
    using AS::nba;

    typedef typename Indirect::Index Index;

  };

  template< typename AS,typename Indirect >
  std::ostream&  operator<< ( std::ostream & os,const SubActiveSet<AS,Indirect>& as );


  /* --- HEAVY CODE --------------------------------------------------------- */
  /* --- HEAVY CODE --------------------------------------------------------- */
  /* --- HEAVY CODE --------------------------------------------------------- */
  template< typename AS,typename Indirect >
  SubActiveSet<AS,Indirect>::
  SubActiveSet( unsigned int nr )
    : AS(nr), self_indirect(1),indirect(self_indirect),isEmpty(true) {}

  template< typename AS,typename Indirect >
  SubActiveSet<AS,Indirect>::
  SubActiveSet( unsigned int nr, Indirect& idx )
    : AS(nr), self_indirect(1),indirect(idx),isEmpty(true) {}

  template< typename AS,typename Indirect >
  SubActiveSet<AS,Indirect>::
  SubActiveSet( const SubActiveSet& clone )
    : AS((const AS&)clone), self_indirect(clone.self_indirect)
    , indirect( clone.ownIndirection()?self_indirect:clone.indirect )
  { }

  template< typename AS,typename Indirect >
  void SubActiveSet<AS,Indirect>::
  reset( void )
  {
    AS::reset();
    isEmpty =true; indirect.resize(0);
  }

  template< typename AS,typename Indirect >
  unsigned int SubActiveSet<AS,Indirect>::
  activeRow( unsigned int ref, Bound::bound_t type )
  {
    assert( (isEmpty&&(nba==0)) || (nba == indirect.size()) );
    unsigned int rowup = AS::activeRow(ref,type);
    return pushIndirectBack(rowup);
  }
  template< typename AS,typename Indirect >
  unsigned int SubActiveSet<AS,Indirect>::
  pushIndirectBack( unsigned int rowup )
  {
    assert( rowup<size() );
    if( isEmpty )
      {
	indirect.resize(1); indirect(0)=rowup;
	isEmpty = false ;
      }
    else
      {
	indirect.conservativeResize( nba );
	indirect(nba-1) = rowup;
      }
    assert( nba>0 );
    return nba-1;
  }
  template< typename AS,typename Indirect >
  void SubActiveSet<AS,Indirect>::
  unactiveRow( unsigned int rowrm )
  {
    assert( nba == indirect.size() );
    assert( rowrm<nba ); // nba>0

    const int internalRowrm = indirect(rowrm);
    if( nba==1 )
      { isEmpty = true ; indirect.resize(0); }
    else
      {
	const Index s = nba-rowrm-1;
	indirect.segment( rowrm,s ) = indirect.tail( s );
	indirect.conservativeResize( nba-1 );
      }
    AS::unactiveRow(internalRowrm);
  }
  template< typename AS,typename Indirect >
  unsigned int SubActiveSet<AS,Indirect>::
  mapInv( unsigned int row ) const
  {
    assert( row<nba );
    AS::mapInv( indirect(row) );
  }

  template< typename AS,typename Indirect >
  unsigned int SubActiveSet<AS,Indirect>::
  map( unsigned int cst ) const
  {
    const int row_ = AS::map(cst);
    for( unsigned int row=0;row<nba;++row )
      { if( indirect(row)==row_ )  return row; }
   }

  template< typename AS,typename Indirect >
  void SubActiveSet<AS,Indirect>::
  setInitialActivation( const AS& as0 )
  {
    assert( &as0!=this );
    reset(); unsigned int row=0;
    for( unsigned int i=0;i<as0.size();++i )
      {
	if( as0.isActive(i) )
	  {
	    AS::active(i,as0.whichBound(i),row++);
	  }
      }
    if( nba>0 )
      {
	isEmpty = false;
	if( nba==1 ) { indirect.resize(1); indirect[0] = 0; }
	else
	  indirect = VectorXi::LinSpaced(nba,0,nba-1);
      }
    assert( row == nba );
  }

  template< typename AS,typename Indirect >
  void SubActiveSet<AS,Indirect>::
  disp( std::ostream& os, bool classic ) const
  {
    if( classic )
      {
	os << " [ ";
	for( unsigned int row=0;row<nba;++row )
	  {
	    const unsigned int cst = mapInv(row);
	    if( whichBound(cst) == Bound::BOUND_INF ) os << "-";
	    else if( whichBound(cst) == Bound::BOUND_SUP ) os << "+";
	    os <<cst << " ";
	  }
	os << " ]";
      }
    else
      {
	AS::disp(os,classic);
	os << "Indirect = " << indirect << std::endl;
      }
  }
  template< typename AS,typename Indirect >
  std::ostream& operator<< ( std::ostream & os,const SubActiveSet<AS,Indirect>& as )
  {    as.disp(os); return os;  }

  /* --- Protected --- */

  template< typename AS,typename Indirect >
  void SubActiveSet<AS,Indirect>::
  active( unsigned int ref, Bound::bound_t type, unsigned int row )
  {
    AS::active(ref,type,row);
    pushIndirectBack(row);
  }




} // namespace soth

#endif
