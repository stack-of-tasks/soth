#define SOTH_DEBUG
#define SOTH_DEBUG_MODE 45
//SOTH_OPEN_DEBUG;
#include "soth/debug.h"
namespace soth
{
  class stage__INIT
  {
  public:stage__INIT( void ) { }// sotDebugTrace::openFile(); }
  };
  stage__INIT sotSOT_initiator;
};


#include "soth/Stage.hpp"
#include <Eigen/QR>
namespace Eigen
{
  #include "soth/DestructiveColPivQR.h"
}

namespace soth
{

  using std::endl;

  Stage::
  Stage( const MatrixXd & inJ, const bound_vector_t & inbounds,BaseY & inY )
    : J(inJ), bounds(inbounds)
    ,Y(inY)
    ,nr(J.rows()),nc(J.cols())
    ,activeSet(nr),freeML_(nr)
    ,W_(nr,nr),ML_(nr,nc),e_(nr)
    ,M(ML_,false),L(ML_,false),e(e_,false)
      //,isWIdenty(true)
    ,W(W_,false)
    ,Ir(L.getRowIndices()),Irn(M.getRowIndices() )
    ,sizeM(0),sizeL(0)
  {
    assert( bounds.size() == J.rows() );
    std::fill( freeML_.begin(),freeML_.end(),true );
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
	/*TODO*/assert( false&&"TODO" );
      }

    activeSet=initialIr;
    sotDEBUG(5) << "initIr = " << (MATLAB)(VectorXi)activeSet << std::endl;

    ML_.setZero();
    SubMatrix<MatrixXd,RowPermutation> Jact( J,activeSet );
    Block<MatrixXd> ML = ML_.topRows(activeSet.nbActive()); //(ML_,0,0,activeSet.nbActive(),nc );
    sotDEBUG(15) << "Ja = " << (MATLAB)Jact << std::endl;
    ML = Jact;
    sotDEBUG(15) << "Ja = " << (MATLAB)ML_ << std::endl;
    Y.applyThisOnTheLeft( ML );

    for( unsigned int i=0;i<activeSet.nbActive();++i )
      {
	freeML_[i]=false;
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

  /* TODO: the previous rank is not usefull anymore, since it correspond to the
   * actual rank of Y. */
  unsigned int Stage::
  computeInitialCOD( const unsigned int previousRank,
		     const ActiveSet & initialIr,
		     BaseY& Yinit)
  {
    sotDEBUG(5) << "J = " << (MATLAB)J << std::endl;

    /* Compute ML=J(initIr,:)*Y. */
    computeInitalJY(initialIr);

    W_.setIdentity();

    /* Set the size of M and L. L is supposed full rank yet. */
    /* M=submatrix(ML,1:previousRank); L=submatrix(ML,previousRank+1:end); */
    M.setColRange(0,previousRank);    M.setRowRange(0,sizeA());  sizeM=previousRank;
    L.setColRange(previousRank,nc);   L.setRowRange(0,sizeA());  sizeL=sizeA();
    sotDEBUG(15) << "MY = " << (MATLAB)M << std::endl;
    sotDEBUG(15) << "LY = " << (MATLAB)L << std::endl;
    sotDEBUG(5) << "e = " << (MATLAB)e << std::endl;
    sotDEBUG(25) << "sizesAML = [" << sizeA() << ", " << sizeM << ", " << sizeL << "]." << std::endl;

    /* A=L'; mQR=QR(A); */
    if (L.cols() == 0)
    {
      sotDEBUG(5) << "col size of L is null, skip the end of initialization" << std::endl;
      L.setRowIndices(VectorXi());
      sizeL=0;
      W.setRowRange(0,sizeA());
      W.setColRange(0,sizeA());
      e.setRowRange(0,sizeA());
      return previousRank;
    }

    Transpose<Block<MatrixXd> > subL = ML_.topRightCorner(sizeA(), nc-previousRank).transpose();
    Block<MatrixXd> subY = Yinit.getNextHouseholderEssential();
    Eigen::DestructiveColPivQR<Transpose<Block<MatrixXd> >, Block<MatrixXd> >
      mQR(subL,subY, EPSILON);
    const MatrixXd & R = mQR.matrixR();
    sotDEBUG(25) << "mR = " << (MATLAB)R << std::endl;
    sotDEBUG(25) << "mQ = " << (MATLAB)Y.getHouseholderEssential() << std::endl;
    sotDEBUG(7) << "ML_ = " << (MATLAB)ML_ << std::endl;
    //    for( Index i=0;i<L.rows();++i ) rowL0(i).tail(nc-sizeM-i-1).setZero();

    /* L=triu(mQR'); */
    const VectorXi & P = mQR.colsPermutation().indices();
    L.setRowIndices(P);    M.setRowIndices(P);
    sotDEBUG(7) << "L0 = " << (MATLAB)L << std::endl;
    sotDEBUG(7) << "M0 = " << (MATLAB)M << std::endl;

    W.setColIndices(P);

    //TODO: is it necessary? don't think so ...
    //activeSet.permuteRows(P);
    W.setRowIndices(P);
    e.setRowIndices(P);

    sotDEBUG(5) << "W0 = " << (MATLAB)W << std::endl;
    sotDEBUG(5) << "e0 = " << (MATLAB)e << std::endl;


    /* for i=rank:-1:1
     *   if( L(i,i)!= 0 ) break;
     *     nullifyLineDeficient( i );
     * sizeL = mQR.rank();
     */
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
    Yinit.increaseRank(sizeL);

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

	sotDEBUG(5) << "Widx2 = " << (MATLAB)W.getColIndices() << std::endl;
	sotDEBUG(5) << Ir(i) << "  " << Ir(row) << std::endl;

	sotDEBUG(5) << "W = " << (MATLAB)W << std::endl;
	sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
      }

    removeARowFromL( row );
  }


  /* Remove a row, and commit the changes in M and W. */
  void Stage::
  removeACrossFromW( const unsigned int & row, const unsigned int & col  )
  {
    unsigned int wrowdown = W.getRowIndices()(row);
    unsigned int wcoldown = W.getColIndices()(col);
    sotDEBUG(15) << "wrowdown =" << wrowdown << std::endl;
    sotDEBUG(15) << "wcoldown =" << wcoldown << std::endl;

    W.removeRow(row);	W.removeCol(col);
    e.removeRow(row);

    M.removeRow(col);
    if( col>=sizeN() )
      { L.removeRow(col-sizeN()); sizeL--; }

    activeSet.unactiveRow(wrowdown);
    freeML_[wcoldown]=true;

    sotDEBUG(5) << "W = " << (MATLAB)W << std::endl;
    sotDEBUG(5) << "M = " << (MATLAB)M << std::endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
    sotDEBUG(10) << "Wi = " << (MATLAB)W.getColIndices() << std::endl;
    sotDEBUG(10) << "Mi = " << (MATLAB)Irn << std::endl;
  }

  /* Remove a row of L, and commit the changes in M and W. */
  void Stage::
  removeARowFromL( unsigned int row )
  {
    L.removeRow(row);
    M.pushRowFront(M.removeRow(row+sizeN()));
    W.pushColFront(W.removeCol(row+sizeN()));
    sizeL--;
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
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
	    GivensSequence & Ydown )
  {
    //sotDEBUGPRIOR(+45);
    sotDEBUG(5) << " --- DOWNDATE ---------------------------- " << std::endl;
    unsigned int colToRemove = removeInW( position );
    bool rankDef = colToRemove >= sizeN();
    removeACrossFromW(position,colToRemove);

    if( rankDef )
      { // Full-rank line removed and no rank promotion: resorbe Hessenberg and propagate.
	sotDEBUG(5) << "Whss = " << (MATLAB)W << std::endl;
	sotDEBUG(5) << "Mhss = " << (MATLAB)M << std::endl;
	sotDEBUG(5) << "Lhss = " << (MATLAB)L << std::endl;

	// TODO: use the knowledge that the first colToRemove-sizeN()
	// of L are properly shaped and does not need any Hess.
	regularizeHessenberg(Ydown);
	L.removeCol(sizeL);

	sotDEBUG(5) << "W2 = " << (MATLAB)W << std::endl;
	sotDEBUG(5) << "M2 = " << (MATLAB)M << std::endl;
	sotDEBUG(5) << "L2= " << (MATLAB)L << std::endl;

	return false;
      }
    else
      {
	sotDEBUG(5) << "Nothing to do." << std::endl;

	return true;
      }
  }

  // Return true if the rank decrease operated at the current stage.
  bool Stage::propagateDowndate( GivensSequence & Ydown,
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
	// MLi := MLi Ydown
	RowL MLi = rowML(i);
	Ydown.applyThisOnTheLeft(MLi);
      }
    if( decreasePreviousRank ) return true;
    sotDEBUG(5) << "W0 = " << (MATLAB)W << std::endl;
    sotDEBUG(5) << "M0 = " << (MATLAB)M << std::endl;
    sotDEBUG(5) << "L0 = " << (MATLAB)L << std::endl;
    L.pushColFront( M.popColBack() );
    sizeM--;

    /* Check is one of the M's grown. */
    for( Index i=0;i<sizeN();++i )
      {
	if( std::abs(ML_(Irn(i),sizeM)) > EPSILON )
	  {
	    /* Remove all the non-zero compononent of ML(i+1:end,sizeM). */
	    sotDEBUG(5) << "Found a non zero at "<<i << std::endl;
	    Block<MatrixXd> ML(ML_,0,0,nr,sizeM+1);
	    for( Index j=i+1;j<sizeN();++j )
	      {
		if( std::abs(ML_(Irn(j),sizeM))<=EPSILON ) continue;
		Givensd G1;
		G1.makeGivens(ML_(Irn(i),sizeM),ML_(Irn(j),sizeM));
		ML.applyOnTheLeft( Irn(i),Irn(j),G1.transpose());
		W_.applyOnTheRight( Irn(i),Irn(j),G1);
	      }
	    /* Commute the lines in L. */
	    M.permuteRow(i,sizeN()-1);
	    W.permuteCol(i,sizeN()-1);
	    L.pushRowFront(Irn(sizeN()-1)); sizeL++;
	    sotDEBUG(5) << "M = " << (MATLAB)M << std::endl;
	    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
	    sotDEBUG(5) << "W = " << (MATLAB)W << std::endl;
	    sotDEBUG(5) << "sizeL = " << sizeL << std::endl;

	    return true;
	  }
      }

    /* No rank upgrade, resorbe hessenberg. */
    regularizeHessenberg(Ydown);
    L.popColBack();
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
    return false;
  }

  void Stage::regularizeHessenberg( GivensSequence & Ydown )
  {
    /*
     * for i=i0:rank-1
     *   gr = GR( L(Ir(i),i),L(Ir(i),i+1),i,i+1 );
     *   L = L*GR;
     *   Ydown.push_back( gr );
     */
    for( unsigned int i=0;i<sizeL;++i )
      {
	RowML MLi = rowMrL0(i);
	sotDEBUG(25) << "MLb = " << (MATLAB)rowMrL0(i) << std::endl;
	Givens G1;
	G1.makeGivensAndApply(MLi,sizeM+i,sizeM+i+1);
	sotDEBUG(25) << "MLa = " << (MATLAB)rowMrL0(i) << std::endl;

	for( unsigned r=i+1;r<sizeL;++r )
	  {
	    RowML MLr = rowMrL0(r) ;
	    G1.applyThisOnTheLeft( MLr );
	  }
	Ydown.push(G1);
      }
  }


  /* Rotate W so that W is 1 on position,position and L|position is
   * at worst hessenberg. */
  unsigned int Stage::
  removeInW( const  unsigned int row )
  {
    sotDEBUG(5) << "W0 = " << (MATLAB)W << std::endl;
    sotDEBUG(5) << "M0 = " << (MATLAB)M << std::endl;
    sotDEBUG(5) << "L0 = " << (MATLAB)L << std::endl;

    int col = 0;
    while( std::abs(W(row,col))< EPSILON ) col++;

    for( unsigned int i=col+1;i<sizeA();++i )
      {
	if( std::abs(W(row,i)-1)< EPSILON ) break;

	/* Wt(row,col) VS Wt(row,i) */
	Givensd G1;
	G1.makeGivens(W(row,col),W(row,i));

	W_.applyOnTheRight( Irn(col),Irn(i),G1 );

	const int rs = rowSize(i);
	if( rs>0 )
	  {
	    /* Apply on 2 specific lines of ML, so ML_ is OK. */
	    Block<MatrixXd> ML(ML_,0,0,nr,rs);
	    ML.applyOnTheLeft( Irn(col),Irn(i),G1.transpose());
	  }
  }

    sotDEBUG(5) << "W = " << (MATLAB)W << std::endl;
    sotDEBUG(5) << "M = " << (MATLAB)M << std::endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
    sotDEBUG(5) << "colToRemove = " << col << endl;

    return col;
  }

  /* --- UPDATE ------------------------------------------------------------- */
  /* --- UPDATE ------------------------------------------------------------- */
  /* --- UPDATE ------------------------------------------------------------- */

  unsigned int Stage::
  update( const ConstraintRef & cst,GivensSequence & Yup )
  {
    sotDEBUG(5) << " --- UPDATE ---------------------------- " << std::endl;
    /*
     * Inew = Unused.pop();
     * Row JupY = row(Inew);
     * JupU = Jup*Y;
     * double norm2=0; double rankJ=0;
     * for i=n:-1:1
     *   norm2+=JupY(i)^2;
     *   if norm2!=0 rankJ=i; break;
     *
     * Ir << Inew
     * W(Inew,Inew)=1;
     * if rankJ>sizeM+rank
     *   // Rank increase
     *   for i=rankJ:-1:sizeM+rank+1
     *     Gr gr; gr.init( JupY,i,i-1,0 ); prod(JupY,gr);
     *     Yup.push_back( gr );
     *     return false;
     * else
     *   // No rank increase;
     *   nullifyLineDeficient(Inew);
     *   return true;
     */
    Index wrowup = activeSet.activeRow( cst.first,cst.second );
    Index wcolup = -1;
    for( unsigned int i=0;i<nr;++i ) if( freeML_[i] ) { freeML_[i]= false; wcolup = i; break; }
    assert( (wcolup >= 0)&&(wcolup<nr) );
    assert( (wrowup >= 0)&&(wrowup<nr) );
    sotDEBUG(5) << "wr=" << wrowup << " wc=" << wcolup << endl;

    sotDEBUG(5) << "cst=" << cst.first << " bound=" << cst.second << endl;
    e_(wrowup) = bounds[cst.first].getBound(cst.second);
    sotDEBUG(5) << "bound="<<bounds[cst.first].getBound(cst.second) <<  endl;
    RowML JupY = ML_.row(wcolup);
    JupY = J.row(cst.first); Y.applyThisOnTheLeft(JupY);
    double norm2=0; int rankJ=sizeM;
    for( Index i=nc-1;i>=sizeM;--i )
      {
	norm2+=JupY(i)*JupY(i);
	if( norm2>EPSILON )
	  { rankJ=i+1; break; }
      }
    sotDEBUG(5) << "JupY = " << (MATLAB)JupY << endl;
    sotDEBUG(5) << "rankUp = " << rankJ << endl;

    /* TODO: add a value in e. */

    if( rankJ>sizeM+sizeL )
      { // Rank increase
	/* Remove the tail of JuY. */
	for( Index i=rankJ-1;i>sizeM+sizeL;--i )
	  {
	    sotDEBUG(45) << "% Right-resorb " << i << endl;
	    Givens G1;
	    G1.makeGivensAndApply(JupY,i-1,i);
	    Yup.push(G1);
	  }
	addARow(wrowup,wcolup);
	L.pushColBack(sizeM+sizeL-1);
      }
    else
      { // No rank increase: regularize.
	if( rankJ>sizeM )
	  {
	    addARow(wrowup,wcolup);
	    nullifyLineDeficient(sizeL-1,rankJ-sizeM);
	  }
	else
	  {
	    addARow(wrowup,wcolup,true);
	  }
      }

    sotDEBUG(5) << "W = " << (MATLAB)W << endl;
    sotDEBUG(5) << "M = " << (MATLAB)M << endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << endl;

    return rankJ;
  }


  void Stage::
  propagateUpdate( GivensSequence & Ydown,
		   unsigned int decreasePreviousRank )
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

    bool defDone = decreasePreviousRank<=sizeM;
    if(! defDone )
      {
	M.pushColBack( L.popColFront() );
	sizeM++;
	// L.pushColFront( sizeM+sizeL );
	// sizeL++;
      }

    sotDEBUG(5) << "M = " << (MATLAB)M << endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << endl;

    MatrixXd Yex(nc,nc); Yex.setIdentity(); Ydown.applyThisOnTheLeft(Yex);
    sotDEBUG(5) << (MATLAB)Yex << endl;
    for( unsigned int i=0;i<sizeA();++i )
      {
	RowL MLi = rowML(i);
	Ydown.applyThisOnTheLeftReduced(MLi);
      }
    sotDEBUG(5) << "M = " << (MATLAB)M << endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << endl;

    // sizeM already increased, so sM+sL is the last col of the Hessenberg.
    if( sizeM+sizeL<=decreasePreviousRank )
      { // L increased a column.
	L.pushColBack( sizeM+sizeL );
      }
    else if(! defDone )
      { // rank decrease ongoing...
	const int rdef = decreasePreviousRank-sizeM;
	assert( (rdef>0)&&(rdef<sizeL) );
	nullifyLineDeficient( rdef,rdef-1 );
      }
    else
      { // already lost the rank, nothing to do.
      }

    sotDEBUG(5) << "M = " << (MATLAB)M << endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << endl;
 }

  void Stage::
  addARow( const Index & wrowup,const Index & wcolup,bool deficient )
  {
    // TODO: clean this mess.
     if(deficient)
       {
     	M.pushRowFront( wcolup );
     	W.pushColFront( wcolup );
     	W.pushRowBack( wrowup );
	e.pushRowBack( wrowup );
     	// clean W.
     	W_.row( wrowup ) .setZero();
     	W_.col( wcolup ) .setZero();
     	W_(wrowup,wcolup ) = 1.0;
       }
     else
      {
	sotDEBUG(5) << "W0 = " << (MATLAB)W << endl;
	sotDEBUG(45) << "Widx = " << (MATLAB)W.getColIndices() << endl;
	sotDEBUG(45) << "Lidx = " << (MATLAB)Irn << endl;
	L.pushRowBack( wcolup );
	M.pushRowBack( wcolup );
	W.pushColBack( wcolup );
	W.pushRowBack( wrowup );
	e.pushRowBack( wrowup );
	// clean W.
	W_.row( wrowup ) .setZero();//setConstant(3.33); //setZero();
	W_.col( wcolup ) .setZero();//setConstant(3.33); //setZero();
	W_(wrowup,wcolup ) = 1.0;
	sizeL++;
      }

     sotDEBUG(5) << "W = " << (MATLAB)W << endl;
     sotDEBUG(5) << "M = " << (MATLAB)M << endl;
     sotDEBUG(5) << "L = " << (MATLAB)L << endl;
  }

  /* --- SOLVER ------------------------------------------------------------- */
  /* --- SOLVER ------------------------------------------------------------- */
  /* --- SOLVER ------------------------------------------------------------- */

  /* Zu=Linv*(Ui'*ei-Mi*Yu(1:rai_1,1)); */
  void Stage::solve( VectorXd& Ytu )
  {
    if (sizeL==0)
    {
      sotDEBUG(10) << "size of L is 0, skipping solve" << std::endl;
      return;
    }
    sotDEBUG(5) << "e = " << (MATLAB)e << std::endl;

    VectorBlock<VectorXd> Ue = Ytu.segment( sizeM,sizeL );
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
    Ue -= Mr*Ytu.head(sizeM);

    sotDEBUG(5) << "Uem = " << (MATLAB)Ue << std::endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;
    soth::solveInPlaceWithLowerTriangular(L,Ue);
    sotDEBUG(5) << "LiUe = " << (MATLAB)Ue << std::endl;


    // TODO: solve =   Ytu + Linv( Ue-MYtu ).
  }


  //err = Ju-e = W [M L 0] Y^u - e
  VectorXd Stage::computeErr(const VectorXd& Ytu)
  {
    //TODO : manage temporary memory ?
    VectorXd tmp = M*Ytu.head(sizeM);
    tmp.tail(sizeL) += L.triangularView<Lower>()*Ytu.segment(sizeM, sizeL);
    return W*tmp - e;
  }


  VectorXd Stage::computeRo(const VectorXd& Ytu)
  {
    //TODO : manage temporary memory ?
    assert(Ytu.size() == nc);
    VectorXd tmp = W.transpose()*e - M*Ytu.head(sizeM);
    tmp.tail(sizeL) += L.triangularView<Lower>()*Ytu.segment(sizeM, sizeL);
    return M.transpose()*tmp;
  }



  /* --- ACCESSORS ---------------------------------------------------------- */
  /* --- ACCESSORS ---------------------------------------------------------- */
  /* --- ACCESSORS ---------------------------------------------------------- */

  /* Get line <r> of the matrix [ L 0 .. 0 ]. */
  Stage::RowL Stage::rowL0( const Index r )
  {
    return ML_.row(Ir(r)).tail(nc-sizeM);
  }


  /* Get line <r> of the matrix [ Mr L 0 .. 0 ] (- Mr = M(Ir,:) -)*/
  Stage::RowML Stage::rowMrL0( const Index r )
  {
    return ML_.row(Ir(r));
  }

  /* Get line <r> of the matrix [ M [0;L] 0 ], headed to the non zero part. */
  Stage::RowL Stage::rowML( const Index r )
  {
    return ML_.row(Irn(r)).head(rowSize(r));
  }

  unsigned int Stage::rowSize( const Index r )
  { return (r<sizeN())?sizeM:sizeM+r-sizeN()+1; }

  /* --- TEST RECOMPOSE ----------------------------------------------------- */
  /* --- TEST RECOMPOSE ----------------------------------------------------- */
  /* --- TEST RECOMPOSE ----------------------------------------------------- */

  /* WMLY = [ W*M W(:,1:rank)*L zeros(sizeA,nc-sizeM-sizeL) ]*Y' */
  void Stage::
  recompose( MatrixXd& WMLY ) const
  {
    sotDEBUGIN(5);
    if( sizeA()==0 )
      {
	assert( (sizeL==0)&&(M.rows()==0)&&(L.rows()==0)&&(L.cols()==0)&&(W.rows()==0)&&(W.rows()==0) );
	WMLY.resize(0,nc);
	return;
      }
    WMLY.resize(sizeA(),nc); WMLY.setZero();
    WMLY.block(0,0,sizeA(),sizeM) = W*M;

    sotDEBUG(5) << "UL = " << (MATLAB)WMLY.block(0,sizeM,sizeA(),sizeL) << std::endl;
    sotDEBUG(5) << "W = " << (MATLAB)W << std::endl;
    sotDEBUG(5) << "U = " << (MATLAB)W.block(0,sizeN(),sizeA(),sizeL) << std::endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << std::endl;

    sotDEBUG(25) << "n = " << sizeN() <<" a = "<< sizeA()
		<< " r = "<< sizeL << endl;
    sotDEBUG(25) << "w="<< W.rows()<<"x" << W.cols()
		<< " -- l=" << L.rows() << "x" << L.cols() << endl;

    if (sizeL != 0)
    {
      sotDEBUG(5) << "U_L = " << (MATLAB)(MatrixXd)(W.block(0,sizeN(),sizeA(),sizeL)*L) << std::endl;


      WMLY.block(0,sizeM,sizeA(),sizeL) = W.block(0,sizeN(),sizeA(),sizeL)*L;
    }
    sotDEBUGIN(5);
    sotDEBUG(5) << "WML = " << (MATLAB)WMLY << std::endl;

    Y.applyTransposeOnTheLeft(WMLY);
    sotDEBUG(5) << "WMLY = " << (MATLAB)WMLY << std::endl;
  }

  bool Stage::
  testRecomposition( void ) const
  {
    MatrixXd Jrec; recompose(Jrec);
    MatrixXd Ja_;   SubMatrix<MatrixXd,RowPermutation> Ja = Jactive(Ja_);
    sotDEBUG(15) << "Jrec="<<(MATLAB)Jrec << endl;
    sotDEBUG(15) << "Ja="<<(MATLAB)Ja << endl;
    bool res;
    if( sizeA() ) res = ((Jrec-Ja).norm()<=EPSILON);
    else res = ( (Jrec.cols()==Ja.cols())&&(Jrec.rows()==Jrec.rows()) );
    sotDEBUG(5) <<"% J: Recomposition  " << ((res)?"OK.":"wrong.") << std::endl;

    VectorXd ea_;
    bool vres;
    if( sizeA() ) vres = (e-eactive(ea_)).norm()<=EPSILON;
    else vres = (e.size()==eactive(ea_).size() );

    sotDEBUG(15) << "e="<<(MATLAB)e << endl;
    sotDEBUG(15) << "ea="<<(MATLAB)eactive(ea_) << endl;
    sotDEBUG(5) <<"% e: Recomposition  " << ((vres)?"OK.":"wrong.") << std::endl;

    return res&&vres;
  }

   /* Return a sub matrix containing the active rows of J, in the
   * same order as given by W. J_ is a matrix where th full
   * J is stored (workspace). */
  SubMatrix<MatrixXd,RowPermutation> Stage::
  Jactive( MatrixXd& J_ ) const
  {
    J_.resize(nr,nc); J_.setConstant(-1.11111);
    for( unsigned int i=0;i<nr;++i )
      {
	if( activeSet.isActive(i) )
	  {
	    J_.row(activeSet.where(i)) = J.row(i);
	    sotDEBUG(15) << "where(" << i << ") = " << activeSet.where(i) << endl;
	  }
      }

    return SubMatrix<MatrixXd,RowPermutation> (J_,W.getRowIndices());
  }

  /* Return a sub vector containing the active rows of e, in the
   * same order as given by W. */
  SubMatrix<VectorXd,RowPermutation> Stage::
  eactive( VectorXd& e_ ) const
  {
    e_.resize(nr); e_.setConstant(-1.11111);
    for( unsigned int i=0;i<nr;++i )
      {
	if( activeSet.isActive(i) )
	  {
	    e_(activeSet.where(i)) = bounds[i].getBound(activeSet.whichBound(i));
	  }
      }

    return  SubVectorXd(e_,W.getRowIndices());
  }

  void Stage::
  show( std::ostream& os, unsigned int stageRef, bool check ) const
  {
    sotDEBUGIN(5);

    activeSet.disp(os);
    os << "sizeA = " << sizeA() << endl;

    MatrixXd J_(nr,nc); J_.setConstant(-1.11111);
    VectorXd e_(nr); e_.setConstant(-1.11111);
    for( unsigned int i=0;i<nr;++i )
      {
	if( activeSet.isActive(i) )
	  {
	    J_.row(activeSet.where(i)) = J.row(i);
	    e_(activeSet.where(i)) = bounds[i].getBound(activeSet.whichBound(i));
	    sotDEBUG(15) << "where(" << i << ") = " << activeSet.where(i) << endl;
	  }
      }

    SubMatrix<MatrixXd,RowPermutation> Ja(J_,W.getRowIndices());
    SubVectorXd ea(e_,W.getRowIndices());

    sotDEBUG(5) << "Iw1"<<stageRef<<" = " << (MATLAB)W.getRowIndices() << std::endl;
    sotDEBUG(5) << "Iw2"<<stageRef<<" = " << (MATLAB)W.getColIndices() << std::endl;
    sotDEBUG(25) << "Ie"<<stageRef<<" = " << (MATLAB)e.getRowIndices() << std::endl;
    sotDEBUG(25) << "Il"<<stageRef<<" = " << (MATLAB)L.getRowIndices() << std::endl;
    sotDEBUG(25) << "J"<<stageRef<<"_ = " << (MATLAB)J_ << std::endl;
    sotDEBUG(5) << "e"<<stageRef<<" = " << (MATLAB)ea << std::endl;
    sotDEBUG(25) << "ML"<<stageRef<<"_ = " << (MATLAB)ML_ << std::endl;

    os << "a"<<stageRef<<" = " << (MATLAB)(Indirect)activeSet << std::endl;
    os << "J"<<stageRef<<" = " << (MATLAB)Ja << std::endl;
    os << "e"<<stageRef<<" = " << (MATLAB)e << std::endl;
    os << "W"<<stageRef<<" = " << (MATLAB)W << std::endl;
    os << "M"<<stageRef<<" = " << (MATLAB)M << std::endl;
    os << "L"<<stageRef<<" = " << (MATLAB)L << std::endl;

    if( check )
    {
      MatrixXd Jrec; recompose(Jrec);
      if (Jrec.rows()>0)
      {
        sotDEBUG(5) << "Jrec="<<(MATLAB)Jrec << endl;
        if((Jrec-Ja).norm()>1e-6) os << "Jrec"<<stageRef<<" = " << (MATLAB)Jrec << std::endl;
        else os <<"% Recomposition OK. " << std::endl;
        if((e-ea).norm()<=1e-6) sotDEBUG(5) <<"% Recomposition e OK. " << std::endl;
        else os << "% Recomposition e not OK. " << std::endl;
      }
      else
      {
        sotDEBUG(5) << "Jrec="<<(MATLAB)Jrec << endl;
        os <<"% Recomposition OK. " << std::endl;
        sotDEBUG(5) <<"% Recomposition e OK. " << std::endl;
      }
    }
  }

  ActiveSet Stage::_allRows(0);
  double Stage::EPSILON = 1e-6;

}; // namespace soth
