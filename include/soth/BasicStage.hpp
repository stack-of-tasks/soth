#ifndef __SOTH_BASIC_STAGE__
#define __SOTH_BASIC_STAGE__

#include <Eigen/Core>
#include <list>
#include <string>
#include "soth/SubMatrix.hpp"
#include "soth/solvers.hpp"
#include "soth/Algebra.hpp"
#include "soth/BaseY.hpp"
#include "soth/Bound.hpp"
#include "soth/ActiveSet.hpp"
#include "soth/Givens.hpp"
#include "soth/Allocator.hpp"

namespace soth
{


  typedef Matrix<Bound, Dynamic,1> VectorBound;
  std::ostream& operator<< (std::ostream& os, const VectorBound& t);

  /* --- STAGE -------------------------------------------------------------- */
  class BasicStage
  {
  private:
    typedef Map<MatrixXd> MapXd;
    typedef Map<VectorBound> MapBound;

    MapXd Jmap;
    MapBound boundsMap;

  protected:
    typedef MatrixBase<MapXd> MatrixXdRef;
    typedef MatrixBase<MapBound> VectorBoundRef;

    const MatrixXdRef & J;
    const VectorBoundRef & bound;
    const unsigned int nr,nc; // nr=nbCols(J), nc=nbRows(J).
    const BaseY& Y;

  public:

    /* Constructor from references on the constraint matrix and
     * vector - no copy. */
    BasicStage( const MatrixXd & J, const VectorBound & bounds, const BaseY& Y  );
    /* Constructor from size and data maps. */
    BasicStage( const unsigned int nr, const unsigned int nc,
		const double * Jdata, const Bound * bdata, const BaseY& Y );

    /* Reset the data map from references - no copy. */
    void set( const MatrixXd & J, const VectorBound & bounds );
    /* Reset the data map from pointers. */
    void set( const double * Jdata, const Bound * bdata );

    unsigned int nbConstraints( void ) const { return nr; }

  public: /* For debug purpose, could be remove on RELEASE. */
    std::string name;
    MatrixXd getJ();
    VectorBound getBound();

  };

}; // namespace soth


#endif // #ifndef __SOTH_BASIC_STAGE__
