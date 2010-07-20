#include <Eigen/Core>
#include <iostream>

using namespace Eigen;
using namespace std;


template<class MatrixType >
struct CustomNicolasAccessor
{
  typedef typename ei_traits<MatrixType>::Scalar Scalar;
  static const int idx1[];
  static const int idx2[];

  inline
  const Scalar operator()( int& row, int& col ) const
  {
    return mat.coeff(idx1[row], idx2[col]);
  }

  inline
  Scalar& operator()( int& row, int& col )
  {
    return mat.coeffRef(idx1[row], idx2[col]);
  }

  inline
  const Scalar operator()( int& idx ) const
  {
    return mat.coeff(idx1[idx%9], idx2[idx/9]);
  }

  inline
  Scalar& operator()( int& idx )
  {
    return mat.coeffRef(idx1[idx%9], idx2[idx/9]);
  }



  MatrixType& mat;

  CustomNicolasAccessor( MatrixType& m )
  :mat( m )
  {}
};

template<class MatrixType >
const int CustomNicolasAccessor<MatrixType>::idx1[] = { 0,1,2,3,4,5,6,7,8 };
template<class MatrixType >
const int CustomNicolasAccessor<MatrixType>::idx2[] = { 0,1,2,3,4,5,6,7,8 };

namespace Eigen {
  template<class MatrixType>
    struct ei_functor_traits< CustomNicolasAccessor<MatrixType> >
    {
      enum {
        Cost = 10,
        PacketAccess = false,
        IsRepeatable = false
      };
    };
}


void testSubmatrix()
{

  MatrixXf orig(9,9);
  int k=0;
  for( int i=0; i<9; ++i )
  {
    for( int j=0; j<9; ++j ){
      orig(i,j) = static_cast<float>(k++); }
  }
  cout << orig << endl << endl;

  typedef CustomNicolasAccessor<MatrixXf> MatrixIdxOp;
  MatrixIdxOp acc( orig );
  CwiseNullaryOp<MatrixIdxOp, MatrixXf> Midx = MatrixXf::NullaryExpr( 5, 4, acc );

  cout << Midx;
  cout << endl << endl;
  //Midx(1,1)=-5;
  MatrixXf o2(5,5); o2 = MatrixXf::Ones(5,5);
  cout << o2 << endl << endl;
  cout << o2*Midx << endl<<endl;

  //Midx = MatrixXf::Random(5,4);
}


/*
 CwiseNullaryOp< CustomNullaryOp, Derived > NullaryExpr  	(  	Index   	 rows,
		Index  	cols,
		const CustomNullaryOp &  	func	
	)
 */










/**
 * Petit main pour construire le tableau map, il fait l'hypothèse suivante:
 * quand tu écris en matlab Mat([1,3,2...],1:3:9) ta matrice de départ est 9x9,
 * stockée en column major mode (le mode par défault en Eigen) et en fait
 * ton indice 1, fair référence à la ligne 0 ou la colonne 0.
*/

/*
int main()
{
  Matrix<int,5,1> rows;
  rows << 1,3,2,5,8;
  Vector4i cols;
  cols << 1, 3, 6, 9;

  for( int j=0; j<4; ++j )
  {
    for( int i=0; i<5; ++i )
    {
      cout << (cols[j]-1)*9+rows[i]-1 << ", ";
    }
  }
  cout << endl;
}
*/
