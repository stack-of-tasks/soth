#ifndef __SOTH_BASIC_STAGE__
#define __SOTH_BASIC_STAGE__

#include <Eigen/Core>
#include <list>
#include <string>
#include "soth/Algebra.hpp"
#include "soth/Bound.hpp"
#include "soth/api.hpp"

#include <boost/noncopyable.hpp>
#include <boost/version.hpp>

#ifndef WITHOUT_NOTIFIOR

#if BOOST_VERSION > 105300
#include <boost/signals2.hpp>
#define NOTIFIOR_SIGNAL2
#else
#include <boost/signals.hpp>
#endif

#endif

namespace soth {
class BaseY;

/* --- STAGE -------------------------------------------------------------- */
class SOTH_EXPORT BasicStage : boost::noncopyable {
 private:
  VectorBound boundsInternal;

  typedef Map<const MatrixXd> MapXd;
  typedef Map<const VectorBound> MapBound;
  typedef MatrixXd::Index Index;

  MapXd Jmap;
  MapBound boundsMap;

 protected:
  typedef MapXd MatrixXdRef;
  typedef MapBound VectorBoundRef;

  const MatrixXdRef& J;
  const VectorBoundRef& bounds;
  const Index nr, nc;  // nr=nbCols(J), nc=nbRows(J).
  const BaseY& Y;

 public:
  /* Constructor from references on the constraint matrix and
   * vector - no copy. */
  BasicStage(const MatrixXd& J, const VectorBound& bounds, const BaseY& Y);
  /* Constructor from size and data maps. */
  BasicStage(const Index nr, const Index nc, const double* Jdata,
             const Bound* bdata, const BaseY& Y);
  /* Same as upper, with bdata:=&boundsInternal. */
  BasicStage(const Index nr, const Index nc, const double* Jdata,
             const BaseY& Y);

  /* Reset the data map from references - no copy. */
  void set(const MatrixXd& J, const VectorBound& bounds);
  /* Reset the data map from pointers. */
  void set(const double* Jdata, const Bound* bdata);

  Index nbConstraints(void) const { return nr; }

  /* Return the J row of the <cst> constraint (in the global ppol,
   * not only in the active pool. */
  VectorXd getJrow(const Index& cst) const;
  /* Return the bound-values of constraint <cst>.*/
  Bound getBoundRow(const Index& cst) const;

 public: /* For debug purpose, could be remove on RELEASE. */
  std::string name;
  MatrixXd getJ() const;
  VectorBound getBounds() const;
  VectorBound& getBoundsInternal();

 public: /* Notification, could be removed conditionnaly to the lack of
            boost::signal. */
#ifndef WITHOUT_NOTIFIOR
  typedef boost::function<void(std::string, ConstraintRef, std::string)>
      listener_function_t;
#ifdef NOTIFIOR_SIGNAL2
  boost::signals2::signal<void(std::string, ConstraintRef, std::string)>
      notifior;
#else
  boost::signal<void(std::string, ConstraintRef, std::string)> notifior;
#endif
#else
  inline void notifior(int, int, std::string) {}
#endif
};

}  // namespace soth

#endif  // #ifndef __SOTH_BASIC_STAGE__
