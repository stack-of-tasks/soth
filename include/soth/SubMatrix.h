#ifndef EIGEN_SUB_MATRIX_H
#define EIGEN_SUB_MATRIX_H

enum
{
  RowPermutation,
  ColPermutation,
  RowAndColPermutation
};

template<typename MatrixType, int PermutationType> class SubMatrix;

template<typename MatrixType, int PermutationType>
struct ei_traits<SubMatrix<MatrixType, PermutationType> >
  : ei_traits<MatrixType>
{
  typedef typename ei_nested<MatrixType>::type MatrixTypeNested;
  typedef typename ei_unref<MatrixTypeNested>::type _MatrixTypeNested;
  typedef typename MatrixType::StorageKind StorageKind;
  enum {
    RowsAtCompileTime = (PermutationType==ColPermutation) ? MatrixType::RowsAtCompileTime : Dynamic,
    ColsAtCompileTime = (PermutationType==RowPermutation) ? MatrixType::ColsAtCompileTime : Dynamic,
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime,
    Flags = _MatrixTypeNested::Flags & HereditaryBits,
    CoeffReadCost = _MatrixTypeNested::CoeffReadCost //Todo : check that
  };
};


template<typename MatrixType, int PermutationType = RowAndColPermutation> class SubMatrix
  : public MatrixBase<SubMatrix<MatrixType, PermutationType> >
{
public:

  typedef MatrixBase<SubMatrix> Base;

//using Base::Scalar;
//typedef typename Base::Scalar Scalar;
//typedef typename Base::Index Index;
EIGEN_DENSE_PUBLIC_INTERFACE(SubMatrix)
typedef Matrix<typename Base::Index, MatrixType::RowsAtCompileTime,1> RowIndices;
typedef Matrix<typename Base::Index, MatrixType::ColsAtCompileTime,1> ColIndices;

  inline SubMatrix(const MatrixType& matrix,
    const RowIndices& rowIndices, const ColIndices& colIndices)
    : m_matrix(matrix), m_rowIndices(rowIndices), m_colIndices(colIndices)
  {
  }

  inline SubMatrix(const MatrixType& matrix, bool permuted=true)
    : m_matrix(matrix)
    , m_rowIndices(permuted ? RowIndices::LinSpaced(0, matrix.rows()-1, matrix.rows()) : RowIndices())
    , m_colIndices(permuted ? ColIndices::LinSpaced(0, matrix.cols()-1, matrix.cols()) : ColIndices())
  {
  }

  EIGEN_INHERIT_ASSIGNMENT_OPERATORS(SubMatrix)

  inline Index rows() const { return m_rowIndices.size(); }
  inline Index cols() const { return m_colIndices.size(); }

  inline Scalar& coeffRef(Index row, Index col)
  {
    return m_matrix.coeffRef(m_rowIndices[row], m_colIndices[col]);
  }

  inline const Scalar coeff(Index row, Index col) const
  {
    return m_matrix.coeff(m_rowIndices[row], m_colIndices[col]);
  }

  inline void setRowIndices(const RowIndices& rowIndices)
  {
    ei_assert(rowIndices.size() <= m_matrix.rows());
    m_rowIndices = rowIndices;
  }

  inline void setColIndices(const ColIndices& colIndices)
  {
    ei_assert(colIndices.size() <= m_matrix.cols());
    m_colIndices = colIndices;
  }

  inline const RowIndices& getRowIndices() const
  {
    return m_rowIndices;
  }

  inline const ColIndices& getColIndices() const
  {
    return m_colIndices;
  }

  inline void permuteRow(Index i, Index j)
  {
    ei_assert(i>=0 && i<m_rowIndices.size());
    ei_assert(j>=0 && j<m_rowIndices.size());
    Index tmp = m_rowIndices[i];
    m_rowIndices[i] = m_rowIndices[j];
    m_rowIndices[j] = tmp;
  }

  inline void permuteCol(Index i, Index j)
  {
    ei_assert(i>=0 && i<m_colIndices.size());
    ei_assert(j>=0 && j<m_colIndices.size());
    Index tmp = m_colIndices[i];
    m_colIndices[i] = m_colIndices[j];
    m_colIndices[j] = tmp;
  }
inline Index removeRow( Index irm )
{
  assert( irm<m_rowIndices.size() );
  const Index res = m_rowIndices(irm);
  const Index s = m_rowIndices.size()-irm-1;
  m_rowIndices.segment( irm,s ) = m_rowIndices.tail( s );
  m_rowIndices.conservativeResize( m_rowIndices.size()-1 );
  return res;
}
inline Index popRowBack()
{
  assert( m_rowIndices.size()>0 );
  Index res = m_rowIndices(m_rowIndices.size()-1);
  m_rowIndices.conservativeResize( m_rowIndices.size()-1 );
  return res;
}
inline Index popRowFront()
{
  return removeRow(0);
}
inline void pushRowBack( Index iadd )
{
  assert( iadd<m_matrix.rows() );
  // TODO: assert( find(iadd,m_rowIndices) == -1 );
  const Index s = m_rowIndices.size();
  m_rowIndices.conservativeResize( s+1 );
  m_rowIndices(s) = iadd;
}
inline void pushRowFront( Index iadd )
{
  assert( iadd<m_matrix.rows() );
  // TODO: assert( find(iadd,m_rowIndices) == -1 );
  const Index s = m_rowIndices.size();
  m_rowIndices.conservativeResize( s+1 );
  for( Index i=s;i>0;--i ){ m_rowIndices(i)=m_rowIndices(i-1); }
  m_rowIndices(0) = iadd;
}
inline Index removeCol( Index irm )
{
  assert( (irm<m_colIndices.size())&&(irm>=0) );
  const Index res = m_colIndices(irm);
  const Index s = m_colIndices.size()-irm-1;
  m_colIndices.segment( irm,s ) = m_colIndices.tail( s );
  m_colIndices.conservativeResize( m_colIndices.size()-1 );
  return res;
}
inline Index popColBack()
{
  Index res = m_colIndices(m_colIndices.size()-1);
  assert( m_colIndices.size()>0 );
  m_colIndices.conservativeResize( m_colIndices.size()-1 );
  return res;
}
inline Index popColFront()
{
  return removeCol(0);
}
inline void pushColBack( Index iadd )
{
  assert( iadd<m_matrix.cols() );  assert( 0<=iadd );
  // TODO: assert( find(iadd,m_colIndices) == -1 );
  const Index s = m_colIndices.size();
  m_colIndices.conservativeResize( s+1 );
  m_colIndices(s) = iadd;
}
inline void pushColFront( Index iadd )
{
  assert( iadd<m_matrix.cols() );
  // TODO: assert( find(iadd,m_colIndices) == -1 );
  const Index s = m_colIndices.size();
  m_colIndices.conservativeResize( s+1 );
  for( Index i=s;i>0;--i ){ m_colIndices(i)=m_colIndices(i-1); }
  m_colIndices(0) = iadd;
}

inline void setColRange( Index start, Index end )
{
  assert( start<=end ); assert( start>=0 ); assert( end<=m_matrix.cols() );

  const Index size = end-start;
  switch (size)
  {
    case 0: m_colIndices = VectorXi(); return;
    case 1: m_colIndices.resize(1); m_colIndices << start; return;
    default: m_colIndices = VectorXi::LinSpaced(start,end-1,size); return;
  }
}
inline void setRowRange( Index start, Index end )
{
  assert( start<=end ); assert( start>=0 ); assert( end<=m_matrix.rows() );

  const Index size = end-start;
  switch (size)
  {
    case 0: m_rowIndices = VectorXi(); return;
    case 1: m_rowIndices.resize(1); m_rowIndices << start; return;
    default: m_rowIndices = VectorXi::LinSpaced(start,end-1,size); return;
  }
}

protected:
  const MatrixType& m_matrix;
  RowIndices m_rowIndices;
  ColIndices m_colIndices;
};


template<typename MatrixType> class SubMatrix<MatrixType, ColPermutation>
  : public MatrixBase<SubMatrix<MatrixType, ColPermutation> >
{
public:

  typedef MatrixBase<SubMatrix> Base;
  EIGEN_DENSE_PUBLIC_INTERFACE(SubMatrix) 
  typedef Matrix<typename Base::Index, MatrixType::ColsAtCompileTime,1> ColIndices;

  inline SubMatrix(const MatrixType& matrix, const ColIndices& colIndices)
    : m_matrix(matrix), m_colIndices(colIndices)
  {
  }

  inline SubMatrix(const MatrixType& matrix, bool permuted=true)
    : m_matrix(matrix)
    , m_colIndices(permuted ? ColIndices::LinSpaced(0, matrix.cols()-1, matrix.cols()) : ColIndices())
  {
  }

  EIGEN_INHERIT_ASSIGNMENT_OPERATORS(SubMatrix)

  inline Index rows() const { return m_matrix.rows(); }
  inline Index cols() const { return m_colIndices.size(); }

  inline Scalar& coeffRef(Index row, Index col)
  {
    return m_matrix.coeffRef(row, m_colIndices[col]);
  }

  inline const Scalar coeff(Index row, Index col) const
  {
    return m_matrix.coeff(row, m_colIndices[col]);
  }

  inline void setColIndices(const ColIndices& colIndices)
  {
    ei_assert(colIndices.size() <= m_matrix.cols());
    m_colIndices = colIndices;
  }

  inline const ColIndices& getColIndices() const
  {
    return m_colIndices;
  }

  inline void permuteCol(Index i, Index j)
  {
    ei_assert(i>=0 && i<m_colIndices.size());
    ei_assert(j>=0 && j<m_colIndices.size());
    Index tmp = m_colIndices[i];
    m_colIndices[i] = m_colIndices[j];
    m_colIndices[j] = tmp;
  }
inline Index removeCol( Index irm )
{
  assert( irm<m_colIndices.size() );
  const Index res = m_colIndices(irm);
  const Index s = m_colIndices.size()-irm-1;
  m_colIndices.segment( irm,s ) = m_colIndices.tail( s );
  m_colIndices.conservativeResize( m_colIndices.size()-1 );
}
inline Index popColBack()
{
  Index res = m_colIndices(m_colIndices.size()-1);
  assert( m_colIndices.size()>0 );
  m_colIndices.conservativeResize( m_colIndices.size()-1 );
  return res;
}
inline Index popColFront()
{
  return removeCol(0);
}
inline void pushColBack( Index iadd )
{
  assert( iadd<m_matrix.cols() );
  // TODO: assert( find(iadd,m_colIndices) == -1 );
  const Index s = m_colIndices.size();
  m_colIndices.conservativeResize( s+1 );
  m_colIndices(s) = iadd;
}
inline void pushColFront( Index iadd )
{
  assert( iadd<m_matrix.cols() );
  // TODO: assert( find(iadd,m_colIndices) == -1 );
  const Index s = m_colIndices.size();
  m_colIndices.conservativeResize( s+1 );
  for( Index i=s;i>0;--i ){ m_colIndices(i)=m_colIndices(i-1); }
  m_colIndices(0) = iadd;
}
inline void setColRange( Index start, Index end )
{
  assert( start<=end ); assert( start>=0 ); assert( end<=m_matrix.cols() );

  const Index size = end-start;
  switch (size)
  {
    case 0: m_colIndices = VectorXi(); return;
    case 1: m_colIndices.resize(1); m_colIndices << start; return;
    default: m_colIndices = VectorXi::LinSpaced(start,end-1,size); return;
  }
}

protected:
  const MatrixType& m_matrix;
  ColIndices m_colIndices;
};


  template<typename MatrixType> class SubMatrix<MatrixType, RowPermutation>
  : public MatrixBase<SubMatrix<MatrixType, RowPermutation> >
{
public:

  typedef MatrixBase<SubMatrix> Base;
  EIGEN_DENSE_PUBLIC_INTERFACE(SubMatrix)
  typedef Matrix<typename Base::Index, MatrixType::RowsAtCompileTime,1> RowIndices;

  inline SubMatrix(const MatrixType& matrix,
    const RowIndices& rowIndices)
    : m_matrix(matrix), m_rowIndices(rowIndices)
  {
  }

  inline SubMatrix(const MatrixType& matrix, bool permuted=true)
    : m_matrix(matrix)
    , m_rowIndices(permuted ? RowIndices::LinSpaced(0, matrix.rows()-1, matrix.rows()) : RowIndices())
  {
  }

  EIGEN_INHERIT_ASSIGNMENT_OPERATORS(SubMatrix)

  inline Index rows() const { return m_rowIndices.size(); }
  inline Index cols() const { return m_matrix.cols(); }

  inline Scalar& coeffRef(Index row, Index col)
  {
    return m_matrix.coeffRef(m_rowIndices[row],col);
  }

  inline const Scalar coeff(Index row, Index col) const
  {
    return m_matrix.coeff(m_rowIndices[row], col);
  }

  inline void setRowIndices(const RowIndices& rowIndices)
  {
    ei_assert(rowIndices.size() <= m_matrix.rows());
    m_rowIndices = rowIndices;
  }

  inline const RowIndices& getRowIndices() const
  {
    return m_rowIndices;
  }

  inline void permuteRow(Index i, Index j)
  {
    ei_assert(i>=0 && i<m_rowIndices.size());
    ei_assert(j>=0 && j<m_rowIndices.size());
    Index tmp = m_rowIndices[i];
    m_rowIndices[i] = m_rowIndices[j];
    m_rowIndices[j] = tmp;
  }
inline Index removeRow( Index irm )
{
  assert( irm<m_rowIndices.size() );
  const Index res = m_rowIndices(irm);
  const Index s = m_rowIndices.size()-irm-1;
  m_rowIndices.segment( irm,s ) = m_rowIndices.tail( s );
  m_rowIndices.conservativeResize( m_rowIndices.size()-1 );
  return res;
}
inline Index popRowBack()
{
  assert( m_rowIndices.size()>0 );
  Index res = m_rowIndices(m_rowIndices.size()-1);
  m_rowIndices.conservativeResize( m_rowIndices.size()-1 );
  return res;
}
inline Index popRowFront()
{
  return removeRow(0);
}
inline void pushRowBack( Index iadd )
{
  assert( iadd<m_matrix.rows() );
  // TODO: assert( find(iadd,m_rowIndices) == -1 );
  const Index s = m_rowIndices.size();
  m_rowIndices.conservativeResize( s+1 );
  m_rowIndices(s) = iadd;
}
inline void pushRowFront( Index iadd )
{
  assert( iadd<m_matrix.rows() );
  // TODO: assert( find(iadd,m_rowIndices) == -1 );
  const Index s = m_rowIndices.size();
  m_rowIndices.conservativeResize( s+1 );
  for( Index i=s;i>0;--i ){ m_rowIndices(i)=m_rowIndices(i-1); }
  m_rowIndices(0) = iadd;
}
inline void setRowRange( Index start, Index end )
{
  assert( start<=end ); assert( start>=0 ); assert( end<=m_matrix.rows() );

  const Index size = end-start;
  switch (size)
  {
    case 0: m_rowIndices = VectorXi(); return;
    case 1: m_rowIndices.resize(1); m_rowIndices << start; return;
    default: m_rowIndices = VectorXi::LinSpaced(start,end-1,size); return;
  }
}
protected:
  const MatrixType& m_matrix;
  RowIndices m_rowIndices;
};



#endif // EIGEN_PERMUTED_MATRIX_H
