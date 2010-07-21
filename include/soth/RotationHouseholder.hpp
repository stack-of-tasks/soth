#ifndef __SOTH_ROTATION_HOUSEHOLDER__
#define __SOTH_ROTATION_HOUSEHOLDER__

#include "soth/Algebra.h"
#include <iostream>
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

    template< typename VectorBase >
    RotationHouseholder( const VectorBase& v, const double& f )
      : h(v),factor(f),size(v.size())
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBase)
    }

    template< typename Derived,typename VectorBase >
    RotationHouseholder( const MatrixBase<Derived>& qr, const VectorBase& coeff, Index j )
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBase)
      init(qr,coeff,j);
    }

    // Init from column j of QR-issued matrix qr.
    template< typename Derived,typename VectorBase >
    void init( const MatrixBase<Derived>& qr, const VectorBase& coeff, Index j );

    /* --- operator -- */

    // M := M*H.
    template< typename Derived >
    void multiplyMatrixLeft( MatrixBase<Derived> & M ) const;
    // M := H*M.
    template< typename Derived >
    void multiplyMatrixRight( MatrixBase<Derived> & M ) const;
    // v := H*v = v'*H.
    template< typename VectorBase >
    void multiplyVector( VectorBase & v ) const;

  }; // class RotationHouseholder

  class HouseholderSequence
  {
  public:
    template< typename Derived,typename VectorBase >
    HouseholderSequence( const MatrixBase<Derived> & mQR,
			 const VectorBase & coeff,
			 const typename Derived::Index& rank  );
    HouseholderSequence( void  ) {}

    void reserve( unsigned int size ) { hhList.reserve(size); }
    void append( const HouseholderSequence& hh );

    // v := H1*...*Hn*v = v*Hn*...*H1.
    template< typename VectorBase >
    void applyThisOnVector( VectorBase & v ) const;
    // v := Hn*...*H1*v = v*H1*...*Hn.
    template< typename VectorBase >
    void applyTransposeOnVector( VectorBase & v ) const;

    // M := M*H1*...*Hn.
    template< typename Derived >
    void applyThisOnTheLeft( MatrixBase<Derived> & M ) const;
    // M := M*Hn*...*H1.
    template< typename Derived >
    void applyTransposeOnTheLeft( MatrixBase<Derived> & M ) const;
    // M := H1*...*Hn*M.
    template< typename Derived >
    void applyThisOnTheRight( MatrixBase<Derived> & M ) const;
    // M := Hn*...*H1*M.
    template< typename Derived >
    void applyTransposeOnTheRight( MatrixBase<Derived> & M ) const;

  protected:
    typedef std::vector< RotationHouseholder > base_seq_t;
    typedef base_seq_t::iterator iterator;
    typedef base_seq_t::const_iterator const_iterator;
    typedef base_seq_t::reverse_iterator riterator;
    typedef base_seq_t::const_reverse_iterator const_riterator;
  private:
    std::vector< RotationHouseholder > hhList;
  };


  /* --- HEAVY CODE --------------------------------------------------------- */
  /* --- HEAVY CODE --------------------------------------------------------- */
  /* --- HEAVY CODE --------------------------------------------------------- */
  template< typename Derived,typename VectorBase >
  void RotationHouseholder::
  init( const MatrixBase<Derived>& qr, const VectorBase& coeff, Index j )
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBase)
    assert( qr.diagonalSize() == coeff.size() ); //CHECK :  coeff.size might be smaller in case the matrix is not full rank ?
    assert( qr.cols()>j );
    assert( j>=0 );

    Index n = qr.rows();
    size = n-j-1;
    h = qr.col(j).tail(n-j-1);
    factor = coeff(j);
  }

  template< typename VectorBase >
  void RotationHouseholder::
  multiplyVector( VectorBase & v ) const
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBase)
    typedef typename VectorBase::Index Index;
    /* v -= b*g*(g'*v)
     * with g=[ 0...0 1 h ]. */
    const Index nv = v.size();
    const Index & nh = size;

    /* Compute g'*v. */
    const Index vstart = nv-nh;
    double gv = v(vstart-1);
    if(nh>0) gv+= h.dot( v.tail(nh) );

    /* Compute b*g'*v. */
    gv *= factor;

    /* Apply v -= h*bgv. */
    v(vstart-1) -= gv;
    if(nh>0) v.tail(nh) -= gv*h;
  }

  // M := H*M;
  template< typename Derived >
  void RotationHouseholder::
  multiplyMatrixRight( MatrixBase<Derived> & M ) const
  {
    const Index nc = M.cols();
    for( Index i=0;i<nc;++i )
      {
	Eigen::MatrixXd::ColXpr Jcol = M.col(i);
	multiplyVector( Jcol );
      }
  }

  // M := M*H;
  template< typename Derived >
  void RotationHouseholder::
  multiplyMatrixLeft( MatrixBase<Derived> & M ) const
  {
    const Index nr = M.rows();
    for( Index i=0;i<nr;++i )
      {
	Eigen::MatrixXd::RowXpr Jrow = M.row(i);
	multiplyVector( Jrow );
      }
  }

  /* --- HEAVY CODE --------------------------------------------------------- */
  template< typename Derived, typename VectorBase >
  HouseholderSequence::
  HouseholderSequence( const MatrixBase<Derived> & mQR,
		       const VectorBase & coeff,
		       const typename Derived::Index& rank )
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBase);
    assert( mQR.diagonalSize()>=rank );
    hhList.resize(rank);
    for( int i=0;i<rank;++i )
      hhList[i].init(mQR,coeff,i);
  }

  inline void HouseholderSequence::
  append( const HouseholderSequence& hh )
  {
    hhList.insert( hhList.end(),hh.hhList.begin(),hh.hhList.end() );
  }


  // M := M*H1*...*Hn.
  template< typename Derived >
  void HouseholderSequence::
  applyThisOnTheLeft( MatrixBase<Derived> & M ) const
  {
    typedef typename Derived::Index Index;
    const Index nr = M.rows();
    for( Index i=0;i<nr;++i )
      {
	typename Derived::RowXpr Mrow = M.row(i);
	/* Apply on left vector is equivalent to apply transpose on right. */
	applyTransposeOnVector( Mrow );
      }
  }

  // M := M*Hn*...*H1.
  template< typename Derived >
  void HouseholderSequence::
  applyTransposeOnTheLeft( MatrixBase<Derived> & M ) const
  {
    typedef typename Derived::Index Index;
    const Index nr = M.rows();
    for( Index i=0;i<nr;++i )
      {
	typename Derived::RowXpr Mrow = M.row(i);
	/* Apply on left vector is equivalent to apply transpose on right. */
	applyThisOnVector( Mrow );
      }
  }


  // M := H1*...*Hn*M = [ H*M(:,1) ... H*M(:,nc) ].
  template< typename Derived >
  void HouseholderSequence::
  applyThisOnTheRight( MatrixBase<Derived> & M ) const
  {
    typedef typename Derived::Index Index;
    const Index nr = M.cols();
    for( Index i=0;i<nr;++i )
      {
	typename Derived::ColXpr Mcol = M.col(i);
	applyThisOnVector( Mcol );
      }
  }
  // M := Hn*...*H1*M.
  template< typename Derived >
  void HouseholderSequence::
  applyTransposeOnTheRight( MatrixBase<Derived> & M ) const
  {
    typedef typename Derived::Index Index;
    const Index nr = M.cols();
    for( Index i=0;i<nr;++i )
      {
	typename Derived::ColXpr Mcol = M.col(i);
	applyTransposeOnVector( Mcol );
      }
  }

  // v := H1*...*Hn*v = v*Hn*...*H1.
  template< typename VectorBase >
  void HouseholderSequence::
  applyThisOnVector( VectorBase & v ) const
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBase)
    for( const_riterator iter=hhList.rbegin(); iter!=hhList.rend();++iter )
      {
	iter->multiplyVector(v);
      }
  }
  // v := Hn*...*H1*v = v*H1*...*Hn.
  template< typename VectorBase >
  void HouseholderSequence::
  applyTransposeOnVector( VectorBase & v ) const
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorBase)
    for( const_iterator iter=hhList.begin();iter!=hhList.end();++iter )
      {
	iter->multiplyVector(v);
      }
  }




}; // namespace soth


#endif //#ifndef __SOTH_ROTATION_HOUSEHOLDER__

