#include "soth/Stage.hpp"
#include "soth/solvers.h"
#include <Eigen/QR>

namespace soth
{

  Stage::
  Stage( const MatrixXd & inJ, const bound_vector_t & inbounds,BaseY & inY )
    : J(inJ), bounds(inbounds)
    ,Y(inY)
    ,nr(J.rows()),nc(J.cols())
    ,activeSet(0)
    ,W_(nr,nr),ML_(nr,nc),e_(nr)
    ,M(ML_,false),L(ML_,false),e(e_,false)
    ,isWIdenty(true),W(W_,false)
    ,Ir(L.getRowIndices()),Irn(M.getRowIndices() )
    ,sizeM(0),sizeL(0),sizeN(0),sizeA(0)
  {
    assert( bounds.size() == J.rows() );
    activeSet.reserve(nr);
  }


  /* --- INITIALISATION OF THE COD ------------------------------------------ */
  /* --- INITIALISATION OF THE COD ------------------------------------------ */
  /* --- INITIALISATION OF THE COD ------------------------------------------ */

  /* Compute ML=J(initIr,:)*Y. */
  void Stage::
  computeInitalJY( const ActiveSet & initialIr )
  {
    if( isAllRow(initialIr) ) { computeInitalJY_allRows(); return; }
    if( initialIr.size()==0 )
      {
	std::cerr << "(#" << __LINE__ << "): TODO: initial IR empty." << std::endl;
	/*TODO*/throw "TODO";
      }

    sizeA=initialIr.size();
    activeSet.resize(sizeA);
    for( unsigned int i=0;i<sizeA;++i )
      {
	const Index & idx = initialIr[i].first;

	MatrixXd::RowXpr MLrow = ML_.row(i);
	MLrow = J.row(idx);
	Y.applyThisOnTheLeft( MLrow );

	e_(i) = bounds[idx].getBound( initialIr[i].second );
	activeSet[i] = initialIr[i];
      }
    std::cout << "JY = " << (MATLAB)ML_ << std::endl;
  }
  /* Compute ML=J(:,:)*Y. */
  void Stage::
  computeInitalJY_allRows(void)
  {
    sizeA=0;
    for( unsigned int i=0;i<nr;++i )
      {
	if( bounds[i].getType() != Bound::BOUND_TWIN ) continue;
	activeSet.resize(activeSet.size()+1);

	MatrixXd::RowXpr MLrow = ML_.row(i);
	MLrow = J.row(i);
	Y.applyThisOnTheLeft( MLrow );

	e_(i) = bounds[i].getBound( Bound::BOUND_TWIN );
	activeSet[i] = ConstraintRef( i,Bound::BOUND_TWIN );
	sizeA++;
     }
  }

  unsigned int Stage::
  computeInitialCOD( const unsigned int previousRank,
		     const ActiveSet & initialIr )
  {
    std::cout << "J = " << (MATLAB)J << std::endl;

    /* Compute ML=J(initIr,:)*Y. */
    computeInitalJY(initialIr);

    /* Set the size of M and L. L is supposed full rank yet. */
    sizeL=sizeA;
    sizeM=previousRank;
    //std::cout << "sizesAML = [" << sizeA << ", " << sizeM << ", " << sizeL << "]." << std::endl;

    /* M=submatrix(ML,1:previousRank); L=submatrix(ML,previousRank+1:end); */
    M.setColRange(0,previousRank);    M.setRowRange(0,sizeA);
    L.setColRange(previousRank,nc);   L.setRowRange(0,sizeA);
    e.setRowRange(0,sizeL);
    std::cout << "MY = " << (MATLAB)M << std::endl;
    std::cout << "LY = " << (MATLAB)L << std::endl;
    std::cout << "e = " << (MATLAB)e << std::endl;

    /* A=L'; mQR=QR(A); */
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> mQR(nr,nc);
    mQR.compute(L.transpose());
    const MatrixXd & QR = mQR.matrixQR();
    std::cout << "mQR = " << (MATLAB)QR << std::endl;

    /* L=triu(mQR'); */
    const VectorXi & P = mQR.colsPermutation().indices();
    for( MatrixXd::Index i=0;i<QR.diagonalSize();++i )
      {
	rowL0(P(i)).head(i+1) =  QR.col(i).head(i+1);
	rowL0(P(i)).tail(nc-sizeM-i-1).setZero();
      }
    L.setRowIndices(P);    M.setRowIndices(P);
    std::cout << "L0 = " << (MATLAB)L << std::endl;
    std::cout << "M0 = " << (MATLAB)M << std::endl;

    W.setRowRange(0,sizeL); W.setColIndices(Ir);
    W_.setIdentity(); // DEBUG
    std::cout << "W0 = " << (MATLAB)W << std::endl;

    /* for i=rank:-1:1
     *   if( L(i,i)!= 0 ) break;
     *     nullifyLineDeficient( i );
     * sizeL = mQR.rank();
     */
    const Index rank = mQR.rank();
    while( sizeL>rank )
      {
	/* Nullify the last line of L, which is of size rank. */
	nullifyLineDeficient( sizeL-1,rank );
      }
    L.setColRange(sizeM,sizeM+sizeL);
    std::cout << "L = " << (MATLAB)L << std::endl;
    std::cout << "W = " << (MATLAB)W << std::endl;

    /* Y=Y*Yup; */
    HouseholderSequence Yup( mQR.matrixQR(),mQR.hCoeffs(),rank );
    Y.composeOnTheRight(Yup);

    return previousRank+sizeL;
  }

  /* Nullify the line <row> which is suppose to be length <in_r> (row-1 by
   * default) by given-rotations comming from the left.
   * Apply the same transformation on M, L and W. Remove the line from L
   * and reorder the lines of M and W by the way.
   * row is the position of the line, in_r its length.
   */
  void Stage::
  nullifyLineDeficient( const Index row, const Index in_r )
  {
    /*
      Jd = L.row(r);
      foreach i in rank:-1:1
      if( Jd(i)==0 ) continue;
      gr= GR(L(Ir(i),i),Jd(i),i,r );
      L=gr*L;
      W=W*gr';
      Ir >> r;
      In << r;
    */
    const Index r = (in_r<0)?row-1:in_r;
    //std::cout << "r = " << r << " , row = " << row << std::endl;

    if( isWIdenty )
      {
	isWIdenty = false;
	W_.setIdentity();
      }

    //std::cout << "WLinit = " << (MATLAB)(MatrixXd)(W*L) << std::endl;
    for( Index i=r-1;i>=0;--i )
      {
	Givensd G1;
	G1.makeGivens(L(i,i),L(row,i));
	Block<MatrixXd> ML(ML_,0,0,nr,sizeM+r);
	ML.applyOnTheLeft( Ir(i),Ir(row),G1.transpose());
	W_.applyOnTheRight( Ir(i),Ir(row),G1);

	//std::cout << "W = " << (MATLAB)W << std::endl;
	//std::cout << "L = " << (MATLAB)L << std::endl;
	//std::cout << "WL = " << (MATLAB)(MatrixXd)(W*L) << std::endl;
      }

    L.removeRow(row);
    M.pushRowFront(M.removeRow(row+sizeN));
    W.pushColFront(W.removeCol(row+sizeN));
    sizeL--; sizeN++;
  }



  /* --- SOLVER ------------------------------------------------------------- */
  /* --- SOLVER ------------------------------------------------------------- */
  /* --- SOLVER ------------------------------------------------------------- */

  /* Zu=Linv*(Ui'*ei-Mi*Yu(1:rai_1,1)); */
  void Stage::solve( VectorXd& Yu )
  {
    VectorBlock<VectorXd> Ue = Yu.segment( sizeM,sizeL );
    if( isWIdenty )
      {	  Ue = Transpositions<-1,-1>(W.getColIndices())*(e);      }
    else
      {
	SubMatrixXd U( W_,W.getRowIndices(),Ir );
	Ue = U.transpose()*e;
      }

    SubMatrixXd Mr( ML_,Ir,M.getColIndices() );
    Ue -= Mr*Yu.tail(sizeM);


    std::cout << "ue = " << (MATLAB)Ue << std::endl;
    std::cout << "L = " << (MATLAB)L << std::endl;
    soth::solveInPlaceWithLowerTriangular(L,Ue);
    std::cout << "LiUe = " << (MATLAB)Ue << std::endl;

    //L.triangularView<Lower>().solveInPlace(Ue);
  }



  /* --- ACCESSORS ---------------------------------------------------------- */
  /* --- ACCESSORS ---------------------------------------------------------- */
  /* --- ACCESSORS ---------------------------------------------------------- */

  /* Get line <r> of the matrix [ 0 ... 0 ; L 0 .. 0 ]. */
  Stage::RowL Stage::rowL0( const unsigned int r )
  {
    return ML_.row(Ir(r)).tail(nc-sizeM);
  }


  /* --- TEST RECOMPOSE ----------------------------------------------------- */
  /* --- TEST RECOMPOSE ----------------------------------------------------- */
  /* --- TEST RECOMPOSE ----------------------------------------------------- */

  /* WMLY = [ W*M W(:,1:rank)*L zeros(sizeA,nc-sizeM-sizeL) ]*Y' */
  void Stage::
  recompose( MatrixXd& WMLY )
  {
    WMLY.resize(sizeA,nc); WMLY.setZero();
    WMLY.block(0,0,sizeA,sizeM) = W*M;
    WMLY.block(0,sizeM,sizeA,sizeL) = W.block(0,sizeN,sizeA,sizeL)*L;
    std::cout << "WML = " << (MATLAB)WMLY << std::endl;

    Y.applyTransposeOnTheLeft(WMLY);
    std::cout << "WMLY = " << (MATLAB)WMLY << std::endl;
  }


  void Stage::
  show( std::ostream& os, unsigned int stageRef, bool check )
  {
    Indirect idx(sizeA);
    for( unsigned int i=0;i<sizeA;++i )
      {	idx(i) = activeSet[i].first;      }
    os << "J"<<stageRef<<" = " << (MATLAB)SubMatrix<MatrixXd,RowPermutation>(J,idx) << std::endl;
    os << "e"<<stageRef<<" = " << (MATLAB)e << std::endl;

    os << "W"<<stageRef<<" = " << (MATLAB)W << std::endl;
    os << "M"<<stageRef<<" = " << (MATLAB)M << std::endl;
    os << "L"<<stageRef<<" = " << (MATLAB)L << std::endl;

    if( check )
      {
	MatrixXd Jrec; recompose(Jrec);
	if((Jrec-J).norm()>1e-6) os << "Jrec"<<stageRef<<" = " << (MATLAB)Jrec << std::endl;
	else os <<"% Recomposition OK. " << std::endl;
      }
  }

  Stage::ActiveSet Stage::_allRows;

}; // namespace soth
