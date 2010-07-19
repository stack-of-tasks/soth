
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <iostream>

namespace boost
{
  namespace numeric
  {
    namespace ublas
    {

      // Index array class
      template<class A>
	class index_array {
	typedef index_array<A> self_type;
      public:
	typedef A array_type;
	typedef const A const_array_type;
	typedef typename A::size_type size_type;
	typedef typename A::difference_type difference_type;
	typedef typename A::value_type value_type;
	typedef typename A::const_reference const_reference;
	typedef typename A::reference reference;
	typedef typename A::const_pointer const_pointer;
	typedef typename A::pointer pointer;

	// Construction and destruction
	BOOST_UBLAS_INLINE
	  index_array ():
	size_ (), data_ () {}
	explicit BOOST_UBLAS_INLINE
	  index_array (size_type size):
	size_ (size), data_ (size) {}
	BOOST_UBLAS_INLINE
	  index_array (size_type size, const array_type &data):
	size_ (size), data_ (data) {}
	BOOST_UBLAS_INLINE
	  index_array (pointer start, pointer stop):
	size_ (stop - start), data_ (stop - start) {
	  std::copy (start, stop, data_.begin ());
	}

	BOOST_UBLAS_INLINE
	  size_type size () const {
	  return data_.size();
	}
	BOOST_UBLAS_INLINE
	  const_array_type data () const {
	  return data_;
	}
	BOOST_UBLAS_INLINE
	  array_type& data () {
	  return data_;
	}

	// Random Access Container
	BOOST_UBLAS_INLINE
	  size_type max_size () const {
	  return size_;
	}

	BOOST_UBLAS_INLINE
	  bool empty () const {
	  return data_.size () == 0;
	}

	// Element access
	BOOST_UBLAS_INLINE
	  const_reference operator () (size_type i) const {
	  BOOST_UBLAS_CHECK (i < size_, bad_index ());
	  return data_ [i];
	}
	BOOST_UBLAS_INLINE
	  reference operator () (size_type i) {
	  BOOST_UBLAS_CHECK (i < size_, bad_index ());
	  return data_ [i];
	}

	BOOST_UBLAS_INLINE
	  const_reference operator [] (size_type i) const {
	  return (*this) (i);
	}
	BOOST_UBLAS_INLINE
	  reference operator [] (size_type i) {
	  return (*this) (i);
	}
	BOOST_UBLAS_INLINE
	  void resize(size_type size, value_type init)
	{
	  data_.resize(size,init);
	  size_ = data_.size();
	}
	BOOST_UBLAS_INLINE
	  void remove( const value_type& val )
	{
	  const size_type psize = size();
	  size_type i;
	  for( i=0;i<psize;++i )
	    if( (*this) (i) == val ) break;
	  for( i=i+1;i<psize;++i )
	    (*this) (i-1) = (*this) (i);
	  resize(psize-1,0);
	}


	// Composition
	BOOST_UBLAS_INLINE
	  index_array compose (const basic_range<size_type, difference_type> &r) const {
	  BOOST_UBLAS_CHECK (r.start () + r.size () <= size_, bad_size ());
	  array_type data (r.size ());
	  for (size_type i = 0; i < r.size (); ++ i)
	    data [i] = data_ [r.start () + i];
	  return index_array (r.size (), data);
	}
	BOOST_UBLAS_INLINE
	  index_array compose (const basic_slice<size_type, difference_type> &s) const {
	  BOOST_UBLAS_CHECK (s.start () + s.stride () * (s.size () - (s.size () > 0)) <= size (), bad_size ());
	  array_type data (s.size ());
	  for (size_type i = 0; i < s.size (); ++ i)
	    data [i] = data_ [s.start () + s.stride () * i];
	  return index_array (s.size (), data);
	}
	BOOST_UBLAS_INLINE
	  index_array compose (const index_array &ia) const {
	  array_type data (ia.size_);
	  for (size_type i = 0; i < ia.size_; ++ i) {
	    BOOST_UBLAS_CHECK (ia.data_ [i] <= size_, bad_size ());
	    data [i] = data_ [ia.data_ [i]];
	  }
	  return index_array (ia.size_, data);
	}

	// Comparison
	template<class OA>
	  BOOST_UBLAS_INLINE
	  bool operator == (const index_array<OA> &ia) const {
	  if (size_ != ia.size_)
	    return false;
	  for (size_type i = 0; i < BOOST_UBLAS_SAME (size_, ia.size_); ++ i)
	    if (data_ [i] != ia.data_ [i])
	      return false;
	  return true;
	}
	template<class OA>
	  BOOST_UBLAS_INLINE
	  bool operator != (const index_array<OA> &ia) const {
	  return ! (*this == ia);
	}

	// Iterator types
      private:
	// Use a index difference
	typedef difference_type const_subiterator_type;

      public:
	typedef indexed_const_iterator<index_array, std::random_access_iterator_tag>
	  const_iterator;

	BOOST_UBLAS_INLINE
	  const_iterator begin () const {
	  return const_iterator (*this, 0);
	}
	BOOST_UBLAS_INLINE
	  const_iterator end () const {
	  return const_iterator (*this, size_);
	}

	// Reverse iterator
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	BOOST_UBLAS_INLINE
	  const_reverse_iterator rbegin () const {
	  return const_reverse_iterator (end ());
	}
	BOOST_UBLAS_INLINE
	  const_reverse_iterator rend () const {
	  return const_reverse_iterator (begin ());
	}

	BOOST_UBLAS_INLINE
	  index_array preprocess (size_type size) const {
	  if (this != &all_)
	    return *this;
	  index_array ia (size);
	  for (size_type i = 0; i < size; ++ i)
	    ia (i) = i;
	  return ia;
	}
	static
	  BOOST_UBLAS_INLINE
	  const index_array &all () {
	  return all_;
	}

      private:
	size_type size_;
	array_type data_;
	static const index_array all_;
      };


      inline index_array< unbounded_array<std::size_t> >&
	operator<< ( index_array< unbounded_array<std::size_t> >& idx, unsigned int i )
	{
	  idx.resize(idx.size()+1,i);
	  return idx;
	}
      inline index_array< unbounded_array<std::size_t> >&
	operator>> ( index_array< unbounded_array<std::size_t> >& idx, unsigned int i )
	{
	  idx.remove(i);
	  return idx;
	}

      template<class A>
	const index_array<A> index_array<A>::all_;

      // Fwd
      template<class A = unbounded_array<std::size_t> >
	class index_array;

    }}} // namespace boost::numeric::ublas;


/* ----BINDING FORTRAN ------------------------------------------------------ */
/* ----BINDING FORTRAN ------------------------------------------------------ */
/* ----BINDING FORTRAN ------------------------------------------------------ */

#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#define LAPACK_DGEQP3 FORTRAN_ID( dgeqp3 )
#define LAPACK_DGEQPF FORTRAN_ID( dgeqpf )
extern "C" {
  void LAPACK_DGEQP3( const int* m, const int* n, double* a, const int* lda, int * jpvt,
                      double* tau, double* work, const int* lwork, int* info );
  void LAPACK_DGEQPF( const int* m, const int* n, double* a, const int* lda, int * jpvt,
                      double* tau, double* work, int* info );
};

namespace boost
{
  namespace numeric
  {
    namespace bindings
    {
      namespace lapack
      {

	template<typename bnuTemplateMatrix>
	  inline int geqp (bnuTemplateMatrix &A,
			   ::boost::numeric::ublas::vector< int >& jp,
			   ::boost::numeric::ublas::vector< double >& tau)
        {
          int const mF= traits::matrix_size1 (A);
          int const nF= traits::matrix_size2 (A);
          if(( nF==0)||(mF==0)) return 0;
          traits::detail::array<double> work(std::max(1, nF*32));

          assert (nF <= traits::vector_size (tau));
          assert (nF <= traits::vector_size (work));

          double* aF =  traits::matrix_storage (A);
          int const ldaF = traits::leading_dimension (A);
          int * jpvtF = traits::vector_storage(jp);
          double * tauF = traits::vector_storage (tau);
          double * workF = traits::vector_storage (work);
          int const lworkF =  traits::vector_size (work);
          int infoF;

          LAPACK_DGEQP3(&mF, &nF, aF, &ldaF, jpvtF, tauF, workF, &lworkF, &infoF);
          //LAPACK_DGEQPF(&mF, &nF, aF, &ldaF, jpvtF, tauF, workF, &infoF);
          return infoF;
        }

      /* template<typename bubTemplateMatrix> */
      /*   inline int geqp (bubTemplateMatrix &A, */
      /*                    ::boost::numeric::ublas::vector< int >& jp, */
      /* 			 ::boost::numeric::ublas::vector<double>& tau) */
      /*   { */
      /*     int const mF= A.size1(); */
      /*     int const nF= A.size2(); */
      /*     if(( nF==0)||(mF==0)) return 0; */
      /*     ::boost::numeric::ublas::vector<double> work(std::max(1, nF*32)); */

      /*     assert (nF <= (int)tau.size()); */
      /*     assert (nF <= (int)work.size()); */

      /* 	  ::boost::numeric::ublas::matrix<double> tmpA; */
      /* 	  tmpA = A; */
      /*     int const ldaF = std::max(A.size1(),A.size2()); */
      /* 	  //traits::leading_dimension (A); */

      /* 	  double* aF =  tmpA.data().begin(); */
      /*     int * jpvtF = jp.data().begin(); */
      /*     double * tauF = tau.data().begin(); */
      /*     double * workF = work.data().begin(); */

      /*     int const lworkF =  work.size(); */
      /*     int infoF; */

      /*     LAPACK_DGEQP3(&mF, &nF, aF, &ldaF, jpvtF, tauF, workF, &lworkF, &infoF); */
      /*     return infoF; */
      /*   } */

      }}}} // End Namespaces


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

namespace soth
{


  namespace bnu = boost::numeric::ublas;
  typedef bnu::matrix<double> MatrixXd;
  typedef bnu::vector<double> VectorXd;
  typedef bnu::index_array<> Indirect;

  typedef bnu::matrix_indirect< MatrixXd,Indirect > SubMatrixXd;
  //typedef bnu::triangular_adaptator< SubMatrixXd > TriSubMatrixXd;


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

    MATLAB( const bnu::vector<double>& v1 );
    MATLAB( const bnu::index_array<>& v1 );
    MATLAB( const bnu::matrix<double>& m1);

    static bool fullPrec;
    std::string str;

  };


}; // namespace soth
