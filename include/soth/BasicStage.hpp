#ifndef __SOTH_BASIC_STAGE__
#define __SOTH_BASIC_STAGE__

#include <Eigen/Core>
#include <list>
#include <string>
#include "soth/Bound.hpp"
#include "soth/Algebra.hpp"

namespace soth
{

  class BaseY;

  /* --- STAGE -------------------------------------------------------------- */
  class BasicStage
  {
  private:
    typedef Map<MatrixXd> MapXd;
    typedef Map<VectorBound> MapBound;

    MapXd Jmap;
    MapBound boundsMap;

  protected:
    typedef MapXd MatrixXdRef;
    typedef MapBound VectorBoundRef;

    const MatrixXdRef & J;
    const VectorBoundRef & bounds;
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
    VectorBound getBounds();

  };

}; // namespace soth


#endif // #ifndef __SOTH_BASIC_STAGE__
