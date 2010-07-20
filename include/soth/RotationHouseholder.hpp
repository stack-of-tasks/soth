#ifndef __SOTH_ROTATION_HOUSEHOLDER__
#define __SOTH_ROTATION_HOUSEHOLDER__


#include <vector>

namespace soth
{

  /* --- ROT HH ------------------------------------------------------------- */
  /* --- ROT HH ------------------------------------------------------------- */
  /* --- ROT HH ------------------------------------------------------------- */

  class RotationHouseholder
  {
  public:
    typedef VectorXd::Index Index;;

  protected:
    VectorXd h;
    double factor;
    Index size;

  public:
    RotationHouseholder( void ) {}
    RotationHouseholder( Index inSize )
      : h(inSize),factor(0),size(inSize) {}

    template< typename VectorGen >
    RotationHouseholder( const VectorGen& v, const double& f )
      : h(v.size()),factor(f),size(v.size())
    {
      for( unsigned int i=0;i<size;++i ) h(i)=v(i);
    }

    template< typename MatrixGen,typename VectorGen >
    RotationHouseholder( const MatrixGen& qr, const VectorGen& coeff, Index j )
    {      init(qr,coeff,j);    }

    // Init from column j of QR-issued matrix qr.
    template< typename MatrixGen,typename VectorGen >
    void init( const MatrixGen& qr, const VectorGen& coeff, Index j );

    /* --- operator -- */

    // M := M*H;
    template< typename MatrixGen >
    MatrixGen & multiplyMatrixLeft( MatrixGen & M );
    // M := H*M;
    template< typename MatrixGen >
    MatrixGen & multiplyMatrixRight( MatrixGen & M );
    // v := H*v = v'*H
    template< typename VectorGen >
    void multiplyVector( VectorGen & v );

  }; // class RotationHouseholder

  typedef std::vector< RotationHouseholder > householder_sequence_t;


  /* --- HEAVY CODE --------------------------------------------------------- */
  template< typename MatrixGen,typename VectorGen >
  void RotationHouseholder::
  init( const MatrixGen& qr, const VectorGen& coeff, Index j )
  {
    assert( qr.diagonalSize() == coeff.size() );
    assert( qr.cols()>j );
    assert( j>=0 );

    Index n = qr.rows();
    size = n-j-1;
    h = qr.col(j).tail(n-j-1);
    factor = coeff(j);
  }

  template< typename VectorGen >
  void RotationHouseholder::
  multiplyVector( VectorGen & v )
  {
    typedef typename VectorGen::Index Index;
    // v -= b*g*(g'*v)
    // with g=[ 0...0 1 h ];
    const Index nv = v.size();
    const Index & nh = size;

    // Compute g'*v
    const Index vstart = nv-nh;
    double gv = v(vstart-1);
    for( Index i=0;i<nh;++i )
      {
	gv += v(vstart+i)*h(i);
      }

    // Compute b*g'*v
    gv *= factor;

    // Apply v -= h*bgv
    v(vstart-1) -= gv;
    for( Index i=0;i<nh;++i )
      {
	 v(vstart+i) -= gv*h(i);
      }

  }

 

  /* template< typename MatrixGen > */
  /* void initFromQR( RotationHouseholder_list_t& hh, MatrixGen& A ) */
  /* { */
  /*   const unsigned int colA=A.size2(); */
  /*   const unsigned int rowA=A.size1(); */
  /*   hh.resize(colA); */
  /*   RotationHouseholder_list_t::iterator iterHH=hh.begin(); */
  /*   for( unsigned int j=0;j<colA;++j,++iterHH ) */
  /*     { */
  /* 	iterHH->init(A,j); */
  /* 	for( unsigned int i=j+1;i<rowA; ++i ) { A(i,j) = 0; } */
  /*     } */
  /* } */




}; // namespace soth


#endif //#ifndef __SOTH_ROTATION_HOUSEHOLDER__

