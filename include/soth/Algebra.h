#include <Eigen/Core>
#include <iostream>


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

namespace soth
{
  //  USING_PART_OF_NAMESPACE_EIGEN;
  using namespace Eigen;

}; // namespace soth




/* -------------------------------------------------------------------------- */
/* --- DEBUG ---------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

namespace soth
{

  template< typename MatrixGen >
    void randMatrix( MatrixGen &A, unsigned int r, unsigned int c )
    {
      A.resize(r,c);
      for( unsigned int i = 0; i<r;++i )
	for( unsigned int j = 0; j<c;++j )
	  { A(i,j) = ((rand()+0.0)/RAND_MAX*2)-1.; }
    }
  template< typename VectorGen >
    void randVector( VectorGen &A, unsigned int r )
    {
      A.resize(r);
      for( unsigned int i = 0; i<r;++i )
	{ A(i) = ((rand()+0.0)/RAND_MAX*2)-1.; }
    }

  struct MATLAB
  {
    friend std::ostream & operator << (std::ostream & os, const MATLAB & m );

    MATLAB( const VectorXd& v1 );
    MATLAB( const VectorXi& v1 );
    MATLAB( const MatrixXd& m1);

    static bool fullPrec;
    std::string str;

  };


}; // namespace soth
