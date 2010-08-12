#ifndef __SOTH_SUB_MATRIX_H__
#define __SOTH_SUB_MATRIX_H__

enum
  {
    RowPermutation,
    ColPermutation,
    RowAndColPermutation
  };

template <typename MatrixType>
struct ei_compute_lvalue_bit
{
  enum {ret = (MatrixType::Flags & (DirectAccessBit | LvalueBit))?LvalueBit:0};
};

template<typename MatrixType, int PermutationType, bool IsSub>
class SubMatrix;

template<typename MatrixType, int PermutationType, bool IsSub>
struct ei_traits<SubMatrix<MatrixType, PermutationType, IsSub> >
  : ei_traits<MatrixType>
{
  typedef typename ei_nested<MatrixType>::type MatrixTypeNested;
  typedef typename ei_unref<MatrixTypeNested>::type _MatrixTypeNested;
  typedef typename MatrixType::StorageKind StorageKind;
  enum {
    RowsAtCompileTime = (PermutationType==ColPermutation) ? MatrixType::RowsAtCompileTime : Dynamic,
    ColsAtCompileTime = (PermutationType==RowPermutation) ? MatrixType::ColsAtCompileTime : Dynamic,
    MaxRowsAtCompileTime = (IsSub ? MatrixType::MaxRowsAtCompileTime : Dynamic),
    MaxColsAtCompileTime = (IsSub ? MatrixType::MaxColsAtCompileTime : Dynamic),
    Flags = _MatrixTypeNested::Flags & HereditaryBits | ei_compute_lvalue_bit<_MatrixTypeNested>::ret,
    CoeffReadCost = _MatrixTypeNested::CoeffReadCost //Todo : check that
  };
};


//-------------------Helpers--------------------//
template<typename VectorType>
struct ei_range_helper
{
  typedef typename VectorType::Index Index;

  static VectorType generate(Index low, Index high)
  {
    Index size = high-low;
    if (size>1)
      return VectorType::LinSpaced(size, low, high-1);
    else if (size==0)
      return VectorType();
    else
      return VectorType::Constant(1,low);
  }
};

template<typename MatrixType,
	 int PermutationType,
	 int rowsAtCompileTime=MatrixType::RowsAtCompileTime,
	 int colsAtCompileTime=MatrixType::ColsAtCompileTime>
struct ei_submatrix_index_helper{};

template<typename MatrixType, int rowsAtCompileTime, int colsAtCompileTime>
struct ei_submatrix_index_helper<MatrixType, RowPermutation, rowsAtCompileTime, colsAtCompileTime>
{
  typedef typename MatrixType::Index Index;
  static Index index(Index row, Index col) {return row;}
};

template<typename MatrixType, int rowsAtCompileTime, int colsAtCompileTime>
struct ei_submatrix_index_helper<MatrixType, ColPermutation, rowsAtCompileTime, colsAtCompileTime>
{
  typedef typename MatrixType::Index Index;
  static Index index(Index row, Index col) {return col;}
};

template<typename MatrixType, int colsAtCompileTime>
struct ei_submatrix_index_helper<MatrixType, RowAndColPermutation, 1, colsAtCompileTime>
{
  typedef typename MatrixType::Index Index;
  static Index index(Index row, Index col) {return col;}
};

template<typename MatrixType, int rowsAtCompileTime>
struct ei_submatrix_index_helper<MatrixType, RowPermutation, rowsAtCompileTime, 1>
{
  typedef typename MatrixType::Index Index;
  static Index index(Index row, Index col) {return row;}
};

//--------------------Memory---------------------//

template<typename MatrixType>
class VoidContainer
{
protected:
  VoidContainer(MatrixType& m) {}
};

template<typename MatrixType>
class MatrixContainer
{
protected:
public:
  MatrixContainer(MatrixType& m) : m_matrix(m) {}
  typename MatrixType::Nested m_matrix;
};

template<typename MatrixType, int PermutationType>
struct ei_choose_container_impl
{
  typedef MatrixContainer<MatrixType> type;
  //  typedef VoidContainer<MatrixType> type;
};

template<typename MatrixType>
struct ei_choose_container_impl<MatrixType, RowAndColPermutation>
{
  typedef MatrixContainer<MatrixType> type;
};



//--------------------Row Implementations---------------------//

template<typename MatrixType, bool IsSub=true>
class NoRowSelectionImpl  // : public MatrixContainer<MatrixType>
{
public:
  typedef typename MatrixType::Index Index;
  typedef Matrix<Index, Dynamic, 1, 0, (IsSub ? MatrixType::MaxRowsAtCompileTime : Dynamic), 1> RowIndices;

protected:
  NoRowSelectionImpl(MatrixType& m, bool defaultPermutation=false) : m_contain(m) {}
  NoRowSelectionImpl(MatrixType& m, RowIndices& indices) :  m_contain(m) {}
  MatrixContainer<MatrixType> m_contain;

public:
  Index rowIndex(Index i) const {return i;}
  Index rows() const {return m_contain.m_matrix.rows();}
};



/* Vector-storred permutation matrix, with simple accessors and modifiors: push, permute, etc. */
template <typename MatrixType, bool IsSub=true>
class RowSelectionImpl
{
public:
  typedef typename MatrixType::Index Index;
  typedef Matrix<Index, Dynamic, 1, 0, (IsSub ? MatrixType::MaxRowsAtCompileTime : Dynamic), 1> RowIndices;

  RowSelectionImpl(MatrixType& m, bool defaultPermutation = false)
    :rowIndices(defaultPermutation?ei_range_helper<RowIndices>::generate(0, m.rows()):RowIndices())
  {
  }

  RowSelectionImpl(MatrixType& m, const RowIndices indices)
    :rowIndices(indices)
  {//TODO ? : check indices, need to know the matrix size
  }

  Index rowIndex(Index i) const
  {
    ei_assert(i>=0 && i<rowIndices.size());
    return rowIndices[i];
  }

  Index rows() const {return rowIndices.size();}

  const RowIndices& getRowIndices() const
  {
    return rowIndices;
  }

  void setRowIndices(const RowIndices& indices)
  {
    //ei_assert(!IsSub || indices.size() < m_matrix.rows());
    rowIndices = indices;
  }

  void setRowRange(Index start, Index end)
  {
    ei_assert(start>=0 && start<rows());
    ei_assert(end>=start && start<rows());
    rowIndices = ei_range_helper<RowIndices>::generate(start, end);
  }

  void permuteRows(Index i, Index j)
  {
    ei_assert(i>=0 && i<rows());
    ei_assert(j>=0 && j<rows());
    Index tmp = rowIndices[i];
    rowIndices[i] = rowIndices[j];
    rowIndices[j] = tmp;
  }

  void pushRowFront(Index index)
  {
    //ei_assert(index<m_matrix.rows());
    const Index s = rows();
    rowIndices.conservativeResize(s+1);
    for (Index i=s; i>0; --i) {rowIndices[i]=rowIndices[i-1];}
    rowIndices[0] = index;
  }

  void pushRowBack(Index index)
  {
    const Index s = rows();
    rowIndices.conservativeResize(s+1);
    rowIndices[s] = index;
  }

  Index removeRow(Index index)
  {
    assert(index<rows());
    const Index res = rowIndices[index];
    const Index s = rows()-index-1;
    rowIndices.segment(index,s) = rowIndices.tail(s);
    rowIndices.conservativeResize(rows()-1);
    return res;
  }

  Index popRowBack()
  {
    assert(rows()>0);
    const Index res = rowIndices[rows()-1];
    rowIndices.conservativeResize(rows()-1);
    return res;
  }

  Index popRowFront()
  {
    assert(rows()>0);
    return removeRow(0);
  }


private:
  RowIndices rowIndices;
};



template<typename MatrixType, int PermutationType, bool IsSub>
struct ei_choose_row_impl
{
  typedef RowSelectionImpl<MatrixType, IsSub> type;
};

template<typename MatrixType, bool IsSub>
struct ei_choose_row_impl<MatrixType, ColPermutation, IsSub>
{
  typedef NoRowSelectionImpl<MatrixType, IsSub> type;
};


//--------------------Col Implementations---------------------//

template<typename MatrixType, bool IsSub=true>
class NoColSelectionImpl // : public MatrixContainer<MatrixType>
{
public:
  typedef typename MatrixType::Index Index;
  typedef Matrix<Index, Dynamic, 1, 0, (IsSub ? MatrixType::MaxColsAtCompileTime : Dynamic), 1> ColIndices;

protected:
  NoColSelectionImpl(MatrixType& m, bool defaultPermutation = false) : m_contain(m) {}
  NoColSelectionImpl(MatrixType& m, ColIndices& indices) : m_contain(m) {}
  MatrixContainer<MatrixType> m_contain;

public:
  Index colIndex(Index i) const {return i;}
  Index cols() const {return m_contain.m_matrix.cols();}
};

template <typename MatrixType, bool IsSub=true>
class ColSelectionImpl
{
public:
  typedef typename MatrixType::Index Index;
  typedef Matrix<Index, Dynamic, 1, 0, (IsSub ? MatrixType::MaxColsAtCompileTime : Dynamic), 1> ColIndices;

  ColSelectionImpl(MatrixType& m, bool defaultPermutation = false)
    : m_contain(m),owned_colIndices(),colIndices( owned_colIndices )
  {
    if( defaultPermutation ) setColRange(0,m.cols());
  }

  /* No copy of the indices, but a storage of the reference. */
  ColSelectionImpl(MatrixType& m,ColIndices& indices)
    : m_contain(m),owned_colIndices(),colIndices( owned_colIndices )
  {}

  /* --- Kernel implementation --- */
  Index colIndex(Index i) const
  {
    ei_assert( inIdxRange(i) );
    ei_assert( inMRange(colIndices[i]) );
    return colIndices[i];
  }
  Index cols() const {return colIndices.size();}
  const ColIndices& getColIndices() const { return colIndices; }
  bool colIndicesOwned() { return &colIndices==owned_colIndices; }

  /* --- Basic setters --- */
  void setColIndices(const ColIndices& indices)
  {
    //ei_assert(!IsSub || indices.size() < m_matrix.cols());
    colIndices = indices;
  }
  void setColRange(Index start, Index end)
  {
    ei_assert( inMRange(start) ); ei_assert( inMRange(std::max(0,end-1)) );
    ei_assert( start<=end );
    colIndices = ei_range_helper<ColIndices>::generate(start, end);
  }
  void permuteCols(Index i, Index j)
  {
    ei_assert(i>=0 && i<cols());
    ei_assert(j>=0 && j<cols());
    Index tmp = colIndices[i];
    colIndices[i] = colIndices[j];
    colIndices[j] = tmp;
  }
  void pushColFront(Index index)
  {
    ei_assert( inMRange(index) );
    const Index s = cols();
    colIndices.conservativeResize(s+1);
    for (Index i=s-1; i>0; --i) {colIndices[i]=colIndices[i-1];}
    colIndices[0] = index;
  }
  void pushColBack(Index index)
  {
    ei_assert( inMRange(index) );
    const Index s = cols();
    colIndices.conservativeResize(s+1);
    colIndices[s] = index;
  }
  Index removeCol(Index index)
  {
    assert(index<cols());
    const Index res = colIndices[index];
    const Index s = cols()-index-1;
    colIndices.segment(index,s) = colIndices.tail(s);
    colIndices.conservativeResize(cols()-1);
    return res;
  }
  Index popColBack()
  {
    assert(cols()>0);
    const Index res = colIndices[cols()-1];
    colIndices.conservativeResize(cols()-1);
    return res;
  }
  Index popColFront()
  {
    assert(cols()>0);
    return removeCol(0);
  }

private:
  MatrixContainer<MatrixType> m_contain;
  ColIndices owned_colIndices;
  ColIndices & colIndices;
  bool inMRange( Index i ) const { return 0<=i && i<m_contain.m_matrix.cols(); }
  bool inIdxRange( Index i ) const { return 0<=i && i<cols(); }
};


template<typename MatrixType, int PermutationType, bool IsSub>
struct ei_choose_col_impl
{
  typedef ColSelectionImpl<MatrixType, IsSub> type;
};

template<typename MatrixType, bool IsSub>
struct ei_choose_col_impl<MatrixType, RowPermutation, IsSub>
{
  typedef NoColSelectionImpl<MatrixType, IsSub> type;
};


/* --- SUB MATRIX CLASS ----------------------------------------------------- */


template<typename MatrixType, int PermutationType = RowAndColPermutation, bool IsSub=true>
class SubMatrix
  : public MatrixBase<SubMatrix<MatrixType, PermutationType, IsSub> >
  , public ei_choose_container_impl<MatrixType, PermutationType>::type
  , public ei_choose_row_impl<MatrixType, PermutationType, IsSub>::type
  , public ei_choose_col_impl<MatrixType, PermutationType, IsSub>::type
{
public:

  typedef MatrixBase<SubMatrix<MatrixType, PermutationType, IsSub> > Base;
  typedef typename ei_choose_container_impl<MatrixType, PermutationType>::type MemoryBase;
  typedef typename ei_choose_row_impl<MatrixType, PermutationType, IsSub>::type RowBase;
  typedef typename ei_choose_col_impl<MatrixType, PermutationType, IsSub>::type ColBase;
  typedef typename RowBase::RowIndices RowIndices;
  typedef typename ColBase::ColIndices ColIndices;

  EIGEN_DENSE_PUBLIC_INTERFACE(SubMatrix)

  typedef Matrix<Index, Dynamic, 1> Indices;

  using RowBase::rows;
  using ColBase::cols;

  inline SubMatrix(MatrixType& matrix, bool defaultPermutation = false )
    : MemoryBase(matrix)
    , RowBase(matrix, defaultPermutation)
    , ColBase(matrix, defaultPermutation)
  {
  }

  inline SubMatrix(MatrixType& matrix, bool defaultPermutationRows, bool defaultPermutationCols )
    : MemoryBase(matrix)
    , RowBase(matrix, defaultPermutationRows)
    , ColBase(matrix, defaultPermutationCols)
  {
  }

  inline SubMatrix(MatrixType& matrix, bool defaultPermutationRows, Indices& indicesCols)
    : MemoryBase(matrix)
    , RowBase(matrix,defaultPermutationRows)
    , ColBase(matrix,indicesCols)
  {
  }
  inline SubMatrix(MatrixType& matrix, Indices& indicesRows, bool defaultPermutationCols)
    : MemoryBase(matrix)
    , RowBase(matrix,indicesRows)
    , ColBase(matrix,defaultPermutationCols)
  {
  }

  inline SubMatrix(MatrixType& matrix, Indices& indicesRows, Indices& indicesCols)
    : MemoryBase(matrix)
    , RowBase(matrix,indicesRows)
    , ColBase(matrix,indicesCols)
  {
  }



  EIGEN_INHERIT_ASSIGNMENT_OPERATORS(SubMatrix)

  inline Scalar& coeffRef(Index row, Index col)
  {
    return MemoryBase::m_matrix.const_cast_derived().coeffRef(rowIndex(row), colIndex(col));
  }

  inline Scalar& coeffRef(Index index)
  {
    return MemoryBase::m_matrix.const_cast_derived()
      .coeffRef(ei_submatrix_index_helper<MatrixType, PermutationType>::index(rowIndex(index), colIndex(index)));
  }

  inline const CoeffReturnType coeff(Index row, Index col) const
  {
    return MemoryBase::m_matrix.coeff(rowIndex(row), colIndex(col));
  }
};

#endif // EIGEN_PERMUTED_MATRIX_H

