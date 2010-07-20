#ifndef __SOTH_BOUND__
#define __SOTH_BOUND__

namespace soth
{
  class Bound
  {
  public:
    enum bound_t
      {
	BOUND_NONE
	,BOUND_INF
	,BOUND_SUP
	,BOUND_DOUBLE
	,BOUND_TWIN // equality
      };

  protected:
    bound_t type;
    // In case of twin bounds, the value is stored in valInf.
    double valInf,valSup;
    double & valTwin;

  public:
    Bound( const double & val, bound_t type );
    Bound( const double & inValInf, const double & inValSup );

    const double& getBound( bound_t type );
    /* Return the bound that is violated, NONE if bound are OK.
     * In case of twin-bounds, no check is performed, NONE is always returned. */
    bound_t check( const double & val );

  }; // Class Bound


  typedef std::vector< Bound > bound_vector_t;

}; // namespace soth

#endif // #ifndef __SOTH_BOUND__
