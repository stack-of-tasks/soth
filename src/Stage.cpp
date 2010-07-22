#include "soth/Stage.hpp"
#include "soth/solvers.h"
#include <Eigen/QR>
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

  /* --- DOWNDATE ----------------------------------------------------------- */
  /* --- DOWNDATE ----------------------------------------------------------- */
  /* --- DOWNDATE ----------------------------------------------------------- */

  // Return true if the rank re-increase operated at the current stage.
  /*
   *   gr = Neutralize row <position> in W <position>
   *   L=gr'*L
   *   bool res;
   *   if( L(In.last(),0) == 0
   *     // Rank deficience: regularize Hessenberg
   *	 Ydown = regularizeHessenberg
   *     res=false;
   *   else
   *     // No rank dec: quit
   *	 Ir << In.pop();
   *	 res=true;
   *   Ir >> r; In >> r;
   *   Unused << r;
   *   return res;
   */
  bool Stage::
  downdate( const unsigned int position,
	    givensd_sequence_t & Ydown )
  {
    std::cout << " --- DOWNDATE ---------------------------- " << std::endl;
    removeInW( position );
    e.removeRow( position );
    /* TODO: remove the component in activeSet. */

    // b. Three possibles cases: rank-deficient line removed, or full-rank remove
    //  and rank promotion, or full-rank removed and no promotion.
    if( position < sizeN ) // Rank-def line removed.
      {
	std::cout << "Nothing to do." << std::endl;
	W.removeRow(position);	W.removeCol(position);
	M.removeRow(position); // No row to remove in L.
	sizeN--; sizeA--;
	return true;
      }
    else if( (sizeN>0)&&(std::abs(ML_( Irn(sizeN-1),sizeM ))>EPSILON) )
      { // Apparition of a none zero coeff on the first deficient L-row.
	// std::cout << "ML_ = " << (MATLAB)ML_ << std::endl;
	// std::cout << "Irn = " << (MATLAB)Irn << std::endl;

	W.removeRow(position);	W.removeCol(position);
	M.removeRow(position);
	std::cout << "W = " << (MATLAB)W << std::endl;
	std::cout << "M = " << (MATLAB)M << std::endl;

	//std::cout << "Lnt = " << (MATLAB)L << std::endl;
	L.removeRow(position-sizeN);
	L.pushRowFront(Irn(sizeN-1));
	std::cout << "L = " << (MATLAB)L << std::endl;
	sizeN--; sizeA--;
	return true;
      }
    else // Full-rank line removed and no rank promotion: resorbe Hessenberg and propagate.
      {
	W.removeRow(position);	W.removeCol(position);
	M.removeRow(position);
	L.removeRow(position-sizeN); sizeL--; sizeA--;
	//std::cout << "Lhss = " << (MATLAB)L << std::endl;
	regularizeHessenberg(Ydown);
	L.removeCol(sizeL);
	sizeL--; sizeA--;
	std::cout << "L = " << (MATLAB)L << std::endl;
	return false;
     }
  }

  // Return true if the rank decrease operated at the current stage.
  bool Stage::propagateDowndate( givensd_sequence_t & Ydown,
				 bool decreasePreviousRank )
  {
    /*
     * M=M*Ydown;
     * if(! decreasePreviousRank ) return true;
     * L.indices2().push_front( M.indice2().pop_back() );
     *
     * foreach i in In
     *   if( L(i,0) == 0 continue;
     *   Ir << i; In >> i;
     *   return true;
     *
     * Ydown += regularizeHessenberg
     * return false
     */

    EI_FOREACH( i,Irn )
      {
	//TODO rowML(i)*=Ydown;
	if( decreasePreviousRank ) return true;
      }
    L.pushColFront( M.popColBack() );
    sizeM--;

    /* Check is one of the M's grown. */
    for( Index i=0;i<sizeN;++i )
      {
	if( std::abs(ML_(Irn(i),sizeM)) > EPSILON )
	  {
	    /* Remove all the non-zero compononent of ML(i+1:end,sizeM). */
	    Block<MatrixXd> ML(ML_,0,0,nr,sizeM+1);
	    for( Index j=i+1;j<sizeN;++j )
	      {
		if( std::abs(ML_(Irn(j),sizeM))<=EPSILON ) continue;
		Givensd G1;
		G1.makeGivens(ML_(Irn(i),sizeM),ML_(Irn(j),sizeM));
		ML.applyOnTheLeft( Irn(i),Irn(j),G1.transpose());
		W_.applyOnTheRight( Irn(i),Irn(j),G1);
	      }
	    /* Commute the lines in L. */
	    L.pushRowFront(Irn(i)); sizeL++;
	    M.permuteRow(i,sizeN-1);
	    std::cout << "M = " << (MATLAB)L << std::endl;
	    std::cout << "L = " << (MATLAB)L << std::endl;

	    return true;
	  }
      }

    /* No rank upgrade, resorbe hessenberg. */
    regularizeHessenberg(Ydown);
    L.popColBack();
    std::cout << "L = " << (MATLAB)L << std::endl;

  }

  void Stage::regularizeHessenberg( givensd_sequence_t & Ydown )
  {
    for( unsigned int i=0;i<sizeL;++i )
      {
	RowML MLi = rowMrL0(i);
	std::cout << "MLi = " << (MATLAB)rowMrL0(i) << std::endl;
	Givensd G1;
	G1.makeGivens(MLi(sizeM+i),MLi(sizeM+i+1),&MLi(sizeM+i));
	MLi(sizeM+i+1)=0;

	// SubMatrix<MatrixXd,RowPermutation> MrL( ML_,Ir );
	// MrL.applyOnTheLeft( sizeM+i,sizeM+i+1,G1 );
	for( unsigned r=i+1;r<sizeL;++r )
	  {
	    rowMrL0(r).applyOnTheRight( sizeM+i,sizeM+i+1,G1 );
	  }

	// TODO: store in Y.
	// Y.
      }
  }


  /*
   * for i=i0:rank-1
   *   gr = GR( L(Ir(i),i),L(Ir(i),i+1),i,i+1 );
   *   L = L*GR;
   *   Ydown.push_back( gr );
   */
  /* Rotate W so that W is 1 on position,position and L|position is at worst hessenberg. */
  void Stage::removeInW( const  unsigned int position )
  {
    // std::cout << "W0 = " << (MATLAB)W << std::endl;
    // std::cout << "M0 = " << (MATLAB)M << std::endl;
    // std::cout << "L0 = " << (MATLAB)L << std::endl;

    for( unsigned int i=0;i<position;++i )
      {
	if( std::abs(W(position,i))< EPSILON ) continue;

	/* Wt(i,position) VS Wt(i+1,position) */
	Givensd G1;
	G1.makeGivens(W(position,i+1),W(position,i));

	W_.applyOnTheRight( Irn(i+1),Irn(i),G1 );

	const int rs = rowSize(i+1);
	assert(rs>0);
	/* Apply on 2 specific lines of ML, so ML_ is OK. */
	Block<MatrixXd> ML(ML_,0,0,nr,rs);
	ML.applyOnTheLeft( Irn(i+1),Irn(i),G1.transpose());
      }

    for( unsigned int i=position+1;i<sizeA;++i )
      {
	if( std::abs(W(position,i))< EPSILON ) continue;

	/* Wt(i,position) VS Wt(position,position) */
	Givensd G1;
	G1.makeGivens(W(position,position),W(position,i));

	W_.applyOnTheRight( Irn(position),Irn(i),G1 );

	const int rs = rowSize(i);
	assert(rs>0);
	/* Apply on 2 specific lines of ML, so ML_ is OK. */
	Block<MatrixXd> ML(ML_,0,0,nr,rs);
	ML.applyOnTheLeft( Irn(position),Irn(i),G1.transpose());
     }

    std::cout << "W = " << (MATLAB)W << std::endl;
    std::cout << "M = " << (MATLAB)M << std::endl;
    std::cout << "L = " << (MATLAB)L << std::endl;
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

    //std::cout << "ue = " << (MATLAB)Ue << std::endl;
    //std::cout << "L = " << (MATLAB)L << std::endl;
    soth::solveInPlaceWithLowerTriangular(L,Ue);
    //std::cout << "LiUe = " << (MATLAB)Ue << std::endl;
  }



  /* --- ACCESSORS ---------------------------------------------------------- */
  /* --- ACCESSORS ---------------------------------------------------------- */
  /* --- ACCESSORS ---------------------------------------------------------- */

  /* Get line <r> of the matrix [ L 0 .. 0 ]. */
  Stage::RowL Stage::rowL0( const unsigned int r )
  {
    return ML_.row(Ir(r)).tail(nc-sizeM);
  }


  /* Get line <r> of the matrix [ Mr L 0 .. 0 ] (- Mr = M(Ir,:) -)*/
  Stage::RowML Stage::rowMrL0( const unsigned int r )
  {
    return ML_.row(Ir(r));
  }

  /* Get line <r> of the matrix [ M [0;L] 0 ], headed to the non zero part. */
  Stage::RowL Stage::rowML( const unsigned int r )
  {
    return ML_.row(Ir(r)).head(rowSize(r));
  }

  unsigned int Stage::rowSize( const unsigned int r )
  { return (r<sizeN)?sizeM:sizeM+r-sizeN; }

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
  double Stage::EPSILON = 1e-6;

}; // namespace soth
