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
#define EI_FOREACH(a,b)  for( Index a=0;a<b.size();++a )


  /* TODO : remove this two functions and replace them in the test by the rand of Eigen. */
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

    MATLAB( const double& x );
    template< typename Derived >
    MATLAB( const MatrixBase<Derived> & m1 )
    {
      if( m1.rows()==0 ) initMatrixRowNull( m1.cols() );
      else if( m1.cols()==0 ) initMatrixColNull( m1.rows() );
      else if( m1.IsVectorAtCompileTime)
	{
	  if( m1.cols()==1 ) initVector( m1.col(0) );
	  else if( m1.rows()==1 ) initVector( m1.row(0) );
	}
      else initMatrix(m1);
    }

    template< typename Derived >
    void initMatrix( const MatrixBase<Derived> & m1 );
    template< typename VectorGen >
    void initVector( const VectorGen & v );
    inline void initMatrixNull( void );
    inline void initMatrixColNull( unsigned int size );
    inline void initMatrixRowNull( unsigned int size );


    static bool fullPrec;
    std::string str;

  };


}; // namespace soth


// --- HEAVY CODE ---
namespace soth
{
  template< typename VectorGen >
    void MATLAB::initVector( const VectorGen & m1 )
    {
      std::ostringstream os; os << "[ ";
      std::ostringstream ostmp;
      for(int i=0;i<m1.size();++i )
	{
	  if( m1[i]<0 ) ostmp << "-"; else ostmp << " ";
	  if(MATLAB::fullPrec||(std::abs(m1[i])>1e-6)) ostmp <<  std::abs(m1[i]);
	  else { ostmp << "0"; }
	  if( m1.size()!=i+1 )
	    {
	      ostmp << ",";
	      const int size = ostmp.str().length();
	      for( unsigned int i=size;i<8;++i) ostmp<<" ";
	      ostmp << "\t";
	    }
	  os << ostmp.str(); ostmp.str("") ;
	}
      os << "  ]';";
      str = os.str();
    }


  void MATLAB::initMatrixNull( void ) { str = "[];"; }
  void MATLAB::initMatrixColNull( unsigned int size )
  {   std::ostringstream os;  os << "zeros("<<size<<",0);";  str = os.str();  }
  void MATLAB::initMatrixRowNull( unsigned int size )
  {   std::ostringstream os;  os << "zeros(0,"<<size<<");";  str = os.str();  }


  template< typename Derived >
    void MATLAB::initMatrix( const MatrixBase<Derived> & m1 )
    {
      std::ostringstream os;
      if( fullPrec ) { os << "[...\n" << m1 << "];"; str=os.str(); return; }
      os << "...\n[ ";
      std::ostringstream ostmp;
      for(int i=0;i<m1.rows();++i )
	{
	  for(int j=0;j<m1.cols();++j )
	    {
	      if( m1(i,j)<0 ) ostmp <<"-"; else ostmp << " ";
	      if(MATLAB::fullPrec||(std::abs(m1(i,j))>1e-6)) ostmp <<  std::abs(m1(i,j));
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


}; // namespace soth





#endif // #ifndef __SOTH_ALGEBRA__
