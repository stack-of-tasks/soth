#ifndef __SOTH_ALGEBRA__
#define __SOTH_ALGEBRA__

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


    template< typename Derived >
    MATLAB( const MatrixBase<Derived> & m1 )
    {
      std::ostringstream os; os << "[ ";
      std::ostringstream ostmp;
      for( unsigned int i=0;i<m1.rows();++i )
	{
	  for( unsigned int j=0;j<m1.cols();++j )
	    {
	      if( m1(i,j)<0 ) ostmp << "-"; else ostmp << " ";
	      if(MATLAB::fullPrec||fabs(m1(i,j))>1e-6) ostmp <<  fabs(m1(i,j));
	      else { ostmp << "0"; }
	      if( m1.cols()!=j+1 )
		{
		  ostmp << ",";
		  const int size = ostmp.str().length();
		  for( unsigned int i=size;i<8;++i) ostmp<<" ";
		  ostmp << "\t";
		}
	      os << ostmp.str(); ostmp.str("") ;
	    }
	  if( m1.rows()!=i+1 ) { os << " ;" << std::endl<<"  "; }
	  else { os << "  ];"; }
	}
      str = os.str();
    }

    static bool fullPrec;
    std::string str;

  };


}; // namespace soth


#endif // #ifndef __SOTH_ALGEBRA__
