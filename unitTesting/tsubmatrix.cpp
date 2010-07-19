#include <Eigen/Core>
#include <iostream>

using namespace Eigen;
using namespace std;


template< class SCALAR, class ORIGINAL_MATRIX >
struct CustomNicolasAccessor
{
  static const int idx1[];
  static const int idx2[];

  inline
  const SCALAR operator()( int& row, int& col ) const
  {
    int coeff = idx2[col]*9+idx1[row];
    return mat[ coeff ];
  }

  inline
  SCALAR& operator()( int& row, int& col )
  {
    int coeff = idx2[col]*9+idx1[row];
    return mat[ coeff ];
  }

  inline
  const SCALAR operator()( int& idx ) const
  {
    return mat[ idx1[idx] ];
  }

  inline
  SCALAR& operator()( int& idx )
  {
    return mat[ idx1[idx] ];
  }



  ORIGINAL_MATRIX& mat;

  CustomNicolasAccessor( ORIGINAL_MATRIX& m )
  :mat( m )
  {}
};

template< class SCALAR, class ORIGINAL_MATRIX >
const int CustomNicolasAccessor<SCALAR,ORIGINAL_MATRIX>::idx1[] = { 0,1,2,3,4,  };
template< class SCALAR, class ORIGINAL_MATRIX >
const int CustomNicolasAccessor<SCALAR,ORIGINAL_MATRIX>::idx2[] = { 0,1, 2, 3  };

namespace Eigen {
  template<typename SCALAR, class ORIGINAL_MATRIX>
    struct ei_functor_traits< CustomNicolasAccessor<SCALAR,ORIGINAL_MATRIX> >
    {
      enum {
        Cost = 10,
        PacketAccess = false,
        IsRepeatable = false
      };
    };
}


int main()
{

  MatrixXf orig(9,9);
  int k=0;
  for( int i=0; i<9; ++i )
  {
    for( int j=0; j<9; ++j ){
      orig(i,j) = k++; }
  }
  cout << orig << endl << endl;

  CustomNicolasAccessor<int, MatrixXf> acc( orig );
  cout << MatrixXf::NullaryExpr( 5, 4, acc );
  cout << endl << endl;
  (MatrixXf::NullaryExpr( 5, 4, acc ))(1,1)=-5;
  MatrixXf o2(5,5); o2 = MatrixXf::Ones(5,5);
  cout << o2 << endl << endl;
  cout << o2*MatrixXf::NullaryExpr( 5, 4, acc ) << endl<<endl;

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
