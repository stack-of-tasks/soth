#define SOTH_DEBUG
#define SOTH_DEBUG_MODE 45
//SOTH_OPEN_DEBUG;
#include "soth/debug.h"
namespace soth
{
  class stage__INIT
  {
  public:stage__INIT( void ) {  sotDebugTrace::openFile();  }
  };
  stage__INIT sotSOT_initiator;
};


#include "soth/Stage.hpp"
#include "soth/solvers.h"
#include <Eigen/QR>
namespace Eigen
{
  #include "soth/DestructiveColPivQR.h"
}

namespace soth
{

  Stage::
  Stage( const MatrixXd & inJ, const bound_vector_t & inbounds,BaseY & inY )
    : J(inJ), bounds(inbounds)
    ,Y(inY)
    ,nr(J.rows()),nc(J.cols())
    ,activeSet(nr)
    ,W_(nr,nr),ML_(nr,nc),e_(nr)
    ,M(ML_,false),L(ML_,false),e(e_,false)
      //,isWIdenty(true)
    ,W(W_,false)
    ,Ir(L.getRowIndices()),Irn(M.getRowIndices() )
    ,sizeM(0),sizeL(0),sizeN(0)
  {
    assert( bounds.size() == J.rows() );
  }


  /* --- INITIALISATION OF THE COD ------------------------------------------ */
  /* --- INITIALISATION OF THE COD ------------------------------------------ */
  /* --- INITIALISATION OF THE COD ------------------------------------------ */

  /* Compute ML=J(initIr,:)*Y. */
  void Stage::
  computeInitalJY( const ActiveSet & initialIr )
  {
    if( isAllRow(initialIr) ) { computeInitalJY_allRows(); return; }
    if( initialIr.nbActive()==0 )
      {
	sotERROR << "TODO: initial IR empty." << std::endl;
	/*TODO*/throw "TODO";
      }

    activeSet=initialIr;
    SubMatrix<MatrixXd,RowPermutation> Jact( J,activeSet );
    Block<MatrixXd> ML = ML_.topRows(activeSet.nbActive()); //(ML_,0,0,activeSet.nbActive(),nc );
    ML = Jact;
    Y.applyThisOnTheLeft( ML );

    for( unsigned int i=0;i<activeSet.nbActive();++i )
      {
	if( activeSet.isActive(i) )
	    e_( activeSet.where(i) ) = bounds[i].getBound( activeSet.whichBound(i) );
      }

    sotDEBUG(15) << "JY = " << (MATLAB)ML_ << std::endl;
    sotDEBUG(15) << "e = " << (MATLAB)e_ << std::endl;
  }
  /* Compute ML=J(:,:)*Y. */
  void Stage::
  computeInitalJY_allRows(void)
  {
    ActiveSet Ir(nr);
    for( unsigned int i=0;i<nr;++i )
      {
	if( bounds[i].getType() != Bound::BOUND_TWIN ) continue;
	int r = Ir.activeRow( i,Bound::BOUND_TWIN );
	assert( r==Ir.nbActive()-1 );
     }
    computeInitalJY( Ir );
  }

  unsigned int Stage::
  computeInitialCOD( const unsigned int previousRank,
		     const ActiveSet & initialIr )
  {
    sotDEBUG(5) << "J = " << (MATLAB)J << std::endl;

    /* Compute ML=J(initIr,:)*Y. */
    computeInitalJY(initialIr);

    /* Set the size of M and L. L is supposed full rank yet. */
    /* M=submatrix(ML,1:previousRank); L=submatrix(ML,previousRank+1:end); */
    M.setColRange(0,previousRank);    M.setRowRange(0,sizeA());  sizeM=previousRank;
    L.setColRange(previousRank,nc);   L.setRowRange(0,sizeA());  sizeL=sizeA();
    e.setRowRange(0,sizeL);
    sotDEBUG(15) << "MY = " << (MATLAB)M << std::endl;
    sotDEBUG(15) << "LY = " << (MATLAB)L << std::endl;
    sotDEBUG(5) << "e = " << (MATLAB)e << std::endl;
    sotDEBUG(25) << "sizesAML = [" << sizeA() << ", " << sizeM << ", " << sizeL << "]." << std::endl;

    /* A=L'; mQR=QR(A); */
    Transpose<Block<MatrixXd> > subL = ML_.topRightCorner(sizeA(), nc-previousRank).transpose();
    Block<MatrixXd> subY = Y.getNextHouseholderEssential();
    Eigen::DestructiveColPivQR<Transpose<Block<MatrixXd> >, Block<MatrixXd> >
      mQR(subL,subY);
    const MatrixXd & R = mQR.matrixR();
    sotDEBUG(25) << "mR = " << (MATLAB)R << std::endl;
    sotDEBUG(25) << "mQ = " << (MATLAB)Y.getHouseholderEssential() << std::endl;

    /* L=triu(mQR'); */
    const VectorXi & P = mQR.colsPermutation().indices();
    L.setRowIndices(P);    M.setRowIndices(P);
    sotDEBUG(7) << "L0 = " << (MATLAB)L << std::endl;
    sotDEBUG(7) << "M0 = " << (MATLAB)M << std::endl;

    W.setRowRange(0,sizeL); W.setColIndices(Ir);
    W_.setIdentity();
    sotDEBUG(5) << "W0 = " << (MATLAB)W << std::endl;

    /* for i=rank:-1:1
     *   if( L(i,i)!= 0 ) break;
     *     nullifyLineDeficient( i );
     * sizeL = mQR.rank();
     */
    sotDEBUG(5) << "VAL = " << L(2,2) << std::endl;
    sotDEBUG(5) << "VAL = " << L.diagonal() << std::endl;
    const Index rank = mQR.rank();
    while( sizeL>rank )
      {
	/* Nullify the last line of L, which is of size rank. */
	sotDEBUG(5) << "Nullify " << sizeL-1 << " / " << rank << std::endl;
    sotDEBUG(5) << "W = " << (MATLAB)W << std::endl;
	nullifyLineDeficient( sizeL-1,rank );
     sotDEBUG(5) << "W = " << (MATLAB)W << std::endl;
     }
    L.setColRange(sizeM,sizeM+sizeL);
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
    sotDEBUG(5) << "W = " << (MATLAB)W << std::endl;

    //sotDEBUG(5) << "check" << std::endl << W* << std::endl;

    /* Y=Y*Yup; */
    //HouseholderSequence Yup( subY,mQR.hCoeffs(),rank );
    //Y.composeOnTheRight(Yup);
    Y.increaseRank(sizeL);

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
    for( Index i=r-1;i>=0;--i )
      {
	if( std::abs(L(row,i))<EPSILON ) continue;
	Givensd G1;
	G1.makeGivens(L(i,i),L(row,i));
	Block<MatrixXd> ML(ML_,0,0,nr,sizeM+r);
	ML.applyOnTheLeft( Ir(i),Ir(row),G1.transpose());
	W_.applyOnTheRight( Ir(i),Ir(row),G1);

	sotDEBUG(5) << "W = " << (MATLAB)W << std::endl;
	sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
	//sotDEBUG(5) << "WL = " << (MATLAB)(MatrixXd)(W*L) << std::endl;
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
    sotDEBUG(5) << " --- DOWNDATE ---------------------------- " << std::endl;
    removeInW( position );
    e.removeRow( position );
    /* TODO: remove the component in activeSet. */

    // b. Three possibles cases: rank-deficient line removed, or full-rank remove
    //  and rank promotion, or full-rank removed and no promotion.
    if( position < sizeN ) // Rank-def line removed.
      {
	sotDEBUG(5) << "Nothing to do." << std::endl;
	W.removeRow(position);	W.removeCol(position);
	M.removeRow(position); // No row to remove in L.
	activeSet.unactiveRow(position);
	sizeN--;
	return true;
      }
    else if( (sizeN>0)&&(std::abs(ML_( Irn(sizeN-1),sizeM ))>EPSILON) )
      { // Apparition of a none zero coeff on the first deficient L-row.
	// sotDEBUG(5) << "ML_ = " << (MATLAB)ML_ << std::endl;
	// sotDEBUG(5) << "Irn = " << (MATLAB)Irn << std::endl;

	W.removeRow(position);	W.removeCol(position);
	M.removeRow(position);
	sotDEBUG(5) << "W = " << (MATLAB)W << std::endl;
	sotDEBUG(5) << "M = " << (MATLAB)M << std::endl;

	//sotDEBUG(5) << "Lnt = " << (MATLAB)L << std::endl;
	L.removeRow(position-sizeN);
	L.pushRowFront(Irn(sizeN-1));
	sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
	activeSet.unactiveRow(position);
	sizeN--;
	return true;
      }
    else // Full-rank line removed and no rank promotion: resorbe Hessenberg and propagate.
      {
	W.removeRow(position);	W.removeCol(position);
	M.removeRow(position);
	L.removeRow(position-sizeN); sizeL--;
	activeSet.unactiveRow(position);
	sotDEBUG(5) << "Lhss = " << (MATLAB)L << std::endl;
	regularizeHessenberg(Ydown);
	L.removeCol(sizeL);
	sizeL--;
	sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
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
      }
    if( decreasePreviousRank ) return true;
    sotDEBUG(5) << "M = " << (MATLAB)M << std::endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
    L.pushColFront( M.popColBack() );
    sizeM--;
    sotDEBUG(5) << "M = " << (MATLAB)M << std::endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;

    /* Check is one of the M's grown. */
    for( Index i=0;i<sizeN;++i )
      {
	if( std::abs(ML_(Irn(i),sizeM)) > EPSILON )
	  {
	    /* Remove all the non-zero compononent of ML(i+1:end,sizeM). */
	    sotDEBUG(5) << "Found a non zero at "<<i << std::endl;
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
	    sotDEBUG(5) << "M = " << (MATLAB)M << std::endl;
	    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;

	    return true;
	  }
      }

    /* No rank upgrade, resorbe hessenberg. */
    regularizeHessenberg(Ydown);
    L.popColBack();
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;

  }

  void Stage::regularizeHessenberg( givensd_sequence_t & Ydown )
  {
    for( unsigned int i=0;i<sizeL;++i )
      {
	RowML MLi = rowMrL0(i);
	sotDEBUG(25) << "MLi = " << (MATLAB)rowMrL0(i) << std::endl;
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
    // sotDEBUG(5) << "W0 = " << (MATLAB)W << std::endl;
    // sotDEBUG(5) << "M0 = " << (MATLAB)M << std::endl;
    // sotDEBUG(5) << "L0 = " << (MATLAB)L << std::endl;

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

    for( unsigned int i=position+1;i<sizeA();++i )
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

    sotDEBUG(5) << "W = " << (MATLAB)W << std::endl;
    sotDEBUG(5) << "M = " << (MATLAB)M << std::endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
  }

  /* --- UPDATE ------------------------------------------------------------- */


  /* --- SOLVER ------------------------------------------------------------- */
  /* --- SOLVER ------------------------------------------------------------- */
  /* --- SOLVER ------------------------------------------------------------- */

  /* Zu=Linv*(Ui'*ei-Mi*Yu(1:rai_1,1)); */
  void Stage::solve( VectorXd& Yu )
  {
    sotDEBUG(5) << "e = " << (MATLAB)e << std::endl;

    VectorBlock<VectorXd> Ue = Yu.segment( sizeM,sizeL );
    // if( isWIdenty )
    //   {	Ue = e; }//  Ue = Transpositions<-1,-1>(W.getColIndices()).transpose()*(e);      }
    // else

    /* TODO: when L0 is full rank, a permuation of e should be enough. */
      {
	SubMatrixXd U( W_,W.getRowIndices(),Ir );
	Ue = U.transpose()*e;
      }

    sotDEBUG(5) << "Ue = " << (MATLAB)Ue << std::endl;
    SubMatrixXd Mr( ML_,Ir,M.getColIndices() );
    Ue -= Mr*Yu.head(sizeM);

    sotDEBUG(5) << "Uem = " << (MATLAB)Ue << std::endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
    soth::solveInPlaceWithLowerTriangular(L,Ue);
    sotDEBUG(5) << "LiUe = " << (MATLAB)Ue << std::endl;
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
    WMLY.resize(sizeA(),nc); WMLY.setZero();
    WMLY.block(0,0,sizeA(),sizeM) = W*M;
    WMLY.block(0,sizeM,sizeA(),sizeL) = W.block(0,sizeN,sizeA(),sizeL)*L;
    sotDEBUG(5) << "WML = " << (MATLAB)WMLY << std::endl;

    Y.applyTransposeOnTheLeft(WMLY);
    sotDEBUG(5) << "WMLY = " << (MATLAB)WMLY << std::endl;
  }


  void Stage::
  show( std::ostream& os, unsigned int stageRef, bool check )
  {
    os << "J"<<stageRef<<" = " << (MATLAB)SubMatrix<MatrixXd,RowPermutation>(J,activeSet) << std::endl;
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

  ActiveSet Stage::_allRows(0);
  double Stage::EPSILON = 1e-6;

}; // namespace soth
