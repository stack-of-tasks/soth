#include "soth/Stage.hpp"
#include <Eigen/QR>

namespace soth
{

  Stage::
  Stage( const MatrixXd & inJ, const bound_vector_t & inbounds,BaseY & inY )
    : J(inJ), bounds(inbounds)
    ,Y(inY)
    ,nr(J.rows()),nc(J.cols())
    ,W_(nr,nr),ML_(nr,nc),e_(nr)
    ,sizeM(0),sizeL(0)
  {
    assert( bounds.size() == J.rows() );
  }


  void Stage::
  computeInitialCOD( const unsigned int previousRank,
		     const Indirect & initialIr )
  {
    std::cout << "J = " << (MATLAB)J << std::endl;

    if( initialIr.size()==0 )
      {
	/*TODO*/throw "TODO";
      }

    /* Compute ML=J(initIr,:)*Y. */
    for( unsigned int i=0;i<initialIr.size();++i )
      {
	MatrixXd::RowXpr MLrow = ML_.row(i);
	MLrow = J.row(initialIr(i));
	Y.applyThisOnTheLeft( MLrow );
      }
    std::cout << "ML = " << (MATLAB)ML_ << std::endl;

    /* Fix the size of M and L. L is supposed full rank yet. */
    sizeA=initialIr.size();
    sizeL=sizeA;
    sizeM=previousRank;
    std::cout << "sizesAML = [" << sizeA << ", " << sizeM << ", "
	      << sizeL << "]." << std::endl;

    // M=submatrix(ML,1:previousRank); L=submatrix(ML,previousRank+1:end);
    Eigen::Block<MatrixXd> M = ML_.block(0,0,sizeA,sizeM);
    std::cout << "M = " << (MATLAB)M << std::endl;
    Eigen::Block<MatrixXd> L = ML_.block(0,sizeM,sizeA,nc-sizeM);
    std::cout << "L = " << (MATLAB)L << std::endl;

    // A=columnMajor(L)  // A==L'
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> mQR(nr,nc);

    // qr(A);
    mQR.compute(L.transpose());
    const MatrixXd & QR = mQR.matrixQR();
    std::cout << "mQR = " << (MATLAB)QR << std::endl;
    for( MatrixXd::Index i=0;i<QR.diagonalSize();++i )
      {
	L.row(i).head(i+1) = QR.col(i).head(i+1);
	L.row(i).tail(nc-sizeM-i-1).setZero();
      }
    std::cout << "L = " << (MATLAB)L << std::endl;

    // RotationHouseHolder_list_t Yup( A );
    HouseholderSequence Yup( mQR.matrixQR(),mQR.hCoeffs() );

    // Y=Y*Yup;
    Y.composeOnTheRight(Yup);

    // for i=rank:-1:1
    //   if( L(i,i)!= 0 ) break;
    // 	nullifyLineDeficient( i );
    sizeL = mQR.rank();

  }




}; // namespace soth
