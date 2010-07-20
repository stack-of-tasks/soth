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
  typedef Matrix<Index, MatrixType::RowsAtCompileTime,1> RowIndices;
  typedef Matrix<Index, MatrixType::ColsAtCompileTime,1> ColIndices;
  EIGEN_DENSE_PUBLIC_INTERFACE(SubMatrix) 

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

  inline const RowIndices& getRowIndices()
  {
    return m_rowIndices;
  }

  inline const ColIndices& getColIndices()
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
  typedef Matrix<Index, MatrixType::ColsAtCompileTime,1> ColIndices;
  EIGEN_DENSE_PUBLIC_INTERFACE(SubMatrix) 

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

  inline const ColIndices& getColIndices()
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

protected:
  const MatrixType& m_matrix;
  ColIndices m_colIndices;
};


  template<typename MatrixType> class SubMatrix<MatrixType, RowPermutation>
  : public MatrixBase<SubMatrix<MatrixType, RowPermutation> >
{
public:

  typedef MatrixBase<SubMatrix> Base;
  typedef Matrix<Index, MatrixType::RowsAtCompileTime,1> RowIndices;
  EIGEN_DENSE_PUBLIC_INTERFACE(SubMatrix)

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

  inline const RowIndices& getRowIndices()
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

protected:
  const MatrixType& m_matrix;
  RowIndices m_rowIndices;
};

#endif // EIGEN_PERMUTED_MATRIX_H
