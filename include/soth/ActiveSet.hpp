#ifndef __SOTH_ASET__
#define __SOTH_ASET__

#include "soth/Algebra.h"
#include "soth/Bound.hpp"

namespace soth
{

  class ActiveSet
  {
  public: /* --- Construction --- */
    ActiveSet( unsigned int nr );
    void reset( void );

  public: /* --- Active set management --- */
    /* Active the constraint <ref> (ie in J(ref,:)) at line <row> (ie in ML(row,:))
     * with bound type. */
    void active( unsigned int ref, Bound::bound_t type, unsigned int row );
    /* Active the given constraint at any free row of ML, and return the
     * position of the selected free row. */
    int activeRow( unsigned int ref, Bound::bound_t type );
    /* Unactive the constraint at line <row> of ML and frees the corresponding line. */
    void unactiveRow( unsigned int row );
    /* Pass the constraint to a twin mode. */
    void freeze( unsigned int ref );


  public: /* --- Accessors --- */
    /* Return the number of active constraint. */
    unsigned int nbActive( void ) const;
    bool isFreezed( unsigned int ref ) const;
    /* Give the reference of the constraint (ie line in J) located at row <row> of
     * the working space ML. */
    unsigned int whichConstraint( unsigned int r ) const;
    unsigned int where( unsigned int ref ) const;
    Bound::bound_t whichBound( unsigned int ref ) const;
    bool isActive( unsigned int ref ) const;
    double sign( unsigned int ref ) const;

  public: /* --- Display --- */
    /* Return a compact of the active line, ordered by row values. */
    void disp( std::ostream& os ) const;

  public: /* --- Deprecated --- */
    /* DEPRECATED*/operator VectorXi (void) const;
    /* DEPRECATED*/void permuteRows( const VectorXi & P );

  protected:
    typedef std::pair<Bound::bound_t,int> cstref_t;
    typedef std::vector<cstref_t> cstref_vector_t;
    cstref_vector_t v;
    std::vector<bool> freerow,freezed;
    unsigned int nba;
  };


} // namespace soth

#endif
