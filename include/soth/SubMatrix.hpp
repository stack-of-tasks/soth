#ifndef __SOTH_SUB_MATRIX_H__
#define __SOTH_SUB_MATRIX_H__

/* The Eigen::MatrixBase derivations have to be done in the
 * Eigen namespace. */
namespace Eigen
{
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
      RowsAtCompileTime = (PermutationType==ColPermutation) ? (MatrixType::RowsAtCompileTime) : Dynamic,
      ColsAtCompileTime = (PermutationType==RowPermutation) ? (MatrixType::ColsAtCompileTime) : Dynamic,
      MaxRowsAtCompileTime = (IsSub ? MatrixType::MaxRowsAtCompileTime : Dynamic),
      MaxColsAtCompileTime = (IsSub ? MatrixType::MaxColsAtCompileTime : Dynamic),
      Flags = (_MatrixTypeNested::Flags & HereditaryBits) | ei_compute_lvalue_bit<_MatrixTypeNested>::ret,
      CoeffReadCost = _MatrixTypeNested::CoeffReadCost //Todo : check that
    };
  };


  /* --- HELPERS ------------------------------------------------------------ */
  /* --- HELPERS ------------------------------------------------------------ */
  /* --- HELPERS ------------------------------------------------------------ */
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

  /* --- CONTAINTER --------------------------------------------------------- */
  /* --- CONTAINTER --------------------------------------------------------- */
  /* --- CONTAINTER --------------------------------------------------------- */

  // TODO: Void containter is not useful anymore, neither is ei_choose_container.
  // Remove them, and propagate the changes in submatrix.
  template<typename MatrixType>
  class VoidContainer
  {
  protected:
    VoidContainer(const MatrixType& m) {}
  };

  template<typename MatrixType>
  class MatrixContainer
  {
  protected:
  public:
    MatrixContainer(const MatrixType& m) : m_matrix(m) {}
    const typename MatrixType::Nested m_matrix;
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


  /* --- ROW INDEX ---------------------------------------------------------- */
  /* --- ROW INDEX ---------------------------------------------------------- */
  /* --- ROW INDEX ---------------------------------------------------------- */

  template<typename MatrixType, bool IsSub=true>
  class NoRowSelectionImpl // : public MatrixContainer<MatrixType>
  {
  public:
    typedef typename MatrixType::Index Index;
    typedef Matrix<Index, Dynamic, 1, 0, (IsSub ? MatrixType::MaxRowsAtCompileTime : Dynamic), 1> RowIndices;

  protected:
    NoRowSelectionImpl(const MatrixType& m, bool defaultPermutation ) : m_contain(m) {}
    NoRowSelectionImpl(const MatrixType& m, const RowIndices indices) : m_contain(m) {}
    NoRowSelectionImpl(const MatrixType& m, RowIndices* indices) : m_contain(m) {}
    MatrixContainer<MatrixType> m_contain;

  public:
    Index rowIndex(Index i) const {return i;}
    Index rows() const {return m_contain.m_matrix.rows();}
  };

  template <typename MatrixType, bool IsSub=true>
  class RowSelectionImpl
  {
  public:
    typedef typename MatrixType::Index Index;
    typedef Matrix<Index, Dynamic, 1, 0, (IsSub ? MatrixType::MaxRowsAtCompileTime : Dynamic), 1> RowIndices;

    RowSelectionImpl(const MatrixType& m, bool defaultPermutation = false)
      : m_contain(m),owned_rowIndices(),rowIndices( owned_rowIndices )
    {
      if( defaultPermutation ) setRowRange(0,m.rows());
    }

    /* No copy of the indices, but a storage of the reference. */
    RowSelectionImpl(const MatrixType& m,const RowIndices indices)
      : m_contain(m),owned_rowIndices(indices),rowIndices(owned_rowIndices )
    { }
    /* No copy of the indices, but a storage of the reference. */
    RowSelectionImpl(const MatrixType& m,RowIndices* indices)
      : m_contain(m),owned_rowIndices(),rowIndices( *indices )
    { }

    /* --- Kernel implementation --- */
    Index rowIndex(Index i) const
    {
      ei_assert( inIdxRange(i) );
      ei_assert( inMRange(rowIndices[i]) );
      return rowIndices[i];
    }
    Index rows() const {return rowIndices.size();}
    const RowIndices& getRowIndices() const { return rowIndices; }
    RowIndices& getRowIndices() { return rowIndices; }
    const Index& getRowIndices(Index i) const { return rowIndices[i]; }
    Index& getRowIndices(Index i) { return rowIndices[i]; }
    bool rowIndicesOwned() { return &rowIndices==owned_rowIndices; }

    /* --- Basic setters --- */
    void setRowIndices(const RowIndices& indices)
    {
      //ei_assert(!IsSub || indices.size() < m_matrix.rows());
      rowIndices = indices;
    }
    void setRowRange(Index start, Index end)
    {
      ei_assert( inMRange(start) ); ei_assert( inMRange(std::max(0,end-1)) );
      ei_assert( start<=end );
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
      ei_assert( inMRange(index) );
      const Index s = rows();
      rowIndices.conservativeResize(s+1);
      /* Cannot use block for this operation. */
      for (Index i=s; i>0; --i) {rowIndices[i]=rowIndices[i-1];}
      rowIndices[0] = index;
    }
    void pushRowBack(Index index)
    {
      ei_assert( inMRange(index) );
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
    MatrixContainer<MatrixType> m_contain;
    RowIndices owned_rowIndices;
    RowIndices & rowIndices;
    bool inMRange( Index i ) const { return 0<=i && i<m_contain.m_matrix.rows(); }
    bool inIdxRange( Index i ) const { return 0<=i && i<rows(); }
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

  /* --- COL INDEX ---------------------------------------------------------- */
  /* --- COL INDEX ---------------------------------------------------------- */
  /* --- COL INDEX ---------------------------------------------------------- */

  template<typename MatrixType, bool IsSub=true>
  class NoColSelectionImpl // : public MatrixContainer<MatrixType>
  {
  public:
    typedef typename MatrixType::Index Index;
    typedef Matrix<Index, Dynamic, 1, 0, (IsSub ? MatrixType::MaxColsAtCompileTime : Dynamic), 1> ColIndices;

  protected:
    NoColSelectionImpl(const MatrixType& m, bool defaultPermutation) : m_contain(m) {}
    NoColSelectionImpl(const MatrixType& m, const ColIndices indices) : m_contain(m) {}
    NoColSelectionImpl(const MatrixType& m, ColIndices* indices) : m_contain(m) {}
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

    ColSelectionImpl(const MatrixType& m, bool defaultPermutation)
      : m_contain(m),owned_colIndices(),colIndices( owned_colIndices )
    {
      if( defaultPermutation ) setColRange(0,m.cols());
    }

    ColSelectionImpl(const MatrixType& m,const ColIndices indices)
      : m_contain(m),owned_colIndices(indices),colIndices( owned_colIndices )
    {}

    /* No copy of the indices, but a storage of the reference. */
    ColSelectionImpl(const MatrixType& m,ColIndices* indices)
      : m_contain(m),owned_colIndices(),colIndices( *indices )
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
    const Index& getColIndices(Index i) const { return colIndices[i]; }
    ColIndices& getColIndices() { return colIndices; }
    Index& getColIndices(Index i) { return colIndices[i]; }
    bool colIndicesOwned() { return &colIndices==owned_colIndices; }

    /* --- Basic setters --- */
    void setColIndices(const ColIndices& indices)
    {
      /* TODO ei_assert extrema of indices. */
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
      for (Index i=s; i>0; --i) {colIndices[i]=colIndices[i-1];}
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
  /* --- SUB MATRIX CLASS ----------------------------------------------------- */
  /* --- SUB MATRIX CLASS ----------------------------------------------------- */
  /* --- SUB MATRIX CLASS ----------------------------------------------------- */

  template <int PermutationType>
  struct ei_choose_assert_selection {};

  template<>
  struct ei_choose_assert_selection<RowPermutation>
  {
    static const int COL_PERMUTATION_IS_NOT_AVAILABLE;
    static const int YOU_SHOULD_HAVE_ONLY_ONE_SUBINDEX;
  };

  template<>
  struct ei_choose_assert_selection<ColPermutation>
  {
    static const int ROW_PERMUTATION_IS_NOT_AVAILABLE;
    static const int YOU_SHOULD_HAVE_ONLY_ONE_SUBINDEX;
  };

  template<>
  struct ei_choose_assert_selection<RowAndColPermutation>
  {
  };


  template<typename MatrixType, int PermutationType = RowAndColPermutation, bool IsSub=true>
  class SubMatrix
    : public MatrixBase<SubMatrix<MatrixType, PermutationType, IsSub> >
    , public ei_choose_container_impl<MatrixType, PermutationType>::type
    , public ei_choose_row_impl<MatrixType, PermutationType, IsSub>::type
    , public ei_choose_col_impl<MatrixType, PermutationType, IsSub>::type
    , public ei_choose_assert_selection<PermutationType>
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
    typedef ei_choose_assert_selection<PermutationType> assert_index;

    inline SubMatrix(MatrixType& matrix )
      : MemoryBase(matrix)
      , RowBase(matrix, false)
      , ColBase(matrix, false)
    {
    }

    inline SubMatrix(MatrixType& matrix, bool defaultPermutation )
      : MemoryBase(matrix)
      , RowBase(matrix, defaultPermutation)
      , ColBase(matrix, defaultPermutation)
    {
      assert(assert_index::YOU_SHOULD_HAVE_ONLY_ONE_SUBINDEX);
    }

    inline SubMatrix(MatrixType& matrix, bool defaultPermutationRows, bool defaultPermutationCols )
      : MemoryBase(matrix)
      , RowBase(matrix, defaultPermutationRows)
      , ColBase(matrix, defaultPermutationCols)
    {
    }

    /* --- By value --- */
    inline SubMatrix(MatrixType& matrix, const Indices indices )
      : MemoryBase(matrix)
      , RowBase(matrix, indices)
      , ColBase(matrix, indices)
    {
      assert(assert_index::YOU_SHOULD_HAVE_ONLY_ONE_SUBINDEX);
    }
    inline SubMatrix(MatrixType& matrix, bool defaultPermutationRows, const Indices indicesCols)
      : MemoryBase(matrix)
      , RowBase(matrix,defaultPermutationRows)
      , ColBase(matrix,indicesCols)
    {
    }
    inline SubMatrix(MatrixType& matrix, const Indices indicesRows, bool defaultPermutationCols)
      : MemoryBase(matrix)
      , RowBase(matrix,indicesRows)
      , ColBase(matrix,defaultPermutationCols)
    {
    }

    inline SubMatrix(MatrixType& matrix, const Indices indicesRows, const Indices indicesCols)
      : MemoryBase(matrix)
      , RowBase(matrix,indicesRows)
      , ColBase(matrix,indicesCols)
    {
    }

    /* --- By reference --- */
    inline SubMatrix(MatrixType& matrix, Indices* indices )
      : MemoryBase(matrix)
      , RowBase(matrix, indices)
      , ColBase(matrix, indices)
    {
      assert(assert_index::YOU_SHOULD_HAVE_ONLY_ONE_SUBINDEX);
    }
    inline SubMatrix(MatrixType& matrix, bool defaultPermutationRows, const Indices* indicesCols)
      : MemoryBase(matrix)
      , RowBase(matrix,defaultPermutationRows)
      , ColBase(matrix,indicesCols)
    {
    }
    inline SubMatrix(MatrixType& matrix, const Indices* indicesRows, bool defaultPermutationCols)
      : MemoryBase(matrix)
      , RowBase(matrix,indicesRows)
      , ColBase(matrix,defaultPermutationCols)
    {
    }

    inline SubMatrix(MatrixType& matrix, Indices* indicesRows, Indices* indicesCols)
      : MemoryBase(matrix)
      , RowBase(matrix,indicesRows)
      , ColBase(matrix,indicesCols)
    {
    }


    /* --- Matrix base interface --- */
    using RowBase::rows;
    using ColBase::cols;

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
    inline const CoeffReturnType coeff(Index index) const
    {
      return MemoryBase::m_matrix.const_cast_derived()
	.coeffRef(ei_submatrix_index_helper<MatrixType, PermutationType>::index(rowIndex(index), colIndex(index)));
    }
  };



  template<typename MatrixType>
  struct ColContainer
  {
    typedef typename MatrixType::ColXpr NestedType;
    NestedType nested;
    typedef typename MatrixType::Index Index;
    ColContainer( MatrixType& m,Index col ) : nested(m.col(col)) {}
  };

  template<typename MatrixType,bool IsSub=true>
  class SubCol
    :protected ColContainer<MatrixType>
    ,public SubMatrix< typename ColContainer<MatrixType>::NestedType,RowPermutation,IsSub >
  {
  public:
    typedef typename MatrixType::Index Index;
    typedef typename ColContainer<MatrixType>::NestedType NestedType;
    typedef SubMatrix< NestedType,RowPermutation,IsSub > SubMatrixBase;
    typedef typename SubMatrixBase::Indices Indices;
    typedef ColContainer<MatrixType> ContainerBase;
    using ContainerBase::nested;

    SubCol( MatrixType& m, const Indices& indicesRows, Index c )
      : ContainerBase(m,c), SubMatrixBase( nested,indicesRows )
    {}

  protected:

  };



} // namespace soth

#endif // __SOTH_SUB_MATRIX_H__

