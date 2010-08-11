#define SOTH_DEBUG
#define SOTH_DEBUG_MODE 45

//SOTH_OPEN_DEBUG;
#include "soth/debug.h"
// namespace soth
// {
//   class stage__INIT
//   {
//   public:stage__INIT( void ) { }// sotDebugTrace::openFile(); }
//   };
//   stage__INIT sotSOT_initiator;
// };


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
    ,W_(nr,nr),ML_(nr,nc),e_(nr),lambda_(nr)
    ,M(ML_,false),L(ML_,false),e(e_,false),lambda(lambda_,false)
      //,isWIdenty(true)
    ,W(W_,false)
    ,Ir(L.getRowIndices()),Irn(M.getRowIndices() )
    ,sizeM(0),sizeL(0)
    ,isReset(false),isInit(false),isOptimumCpt(false),isLagrangeCpt(false)
  {
    assert( bounds.size() == J.rows() );
    std::fill( freeML_.begin(),freeML_.end(),true );
  }


  /* --- INITIALISATION OF THE COD ------------------------------------------ */
  /* --- INITIALISATION OF THE COD ------------------------------------------ */
  /* --- INITIALISATION OF THE COD ------------------------------------------ */

  void Stage::
  reset( void )
  {
    sotDEBUG(45) << "# In {" << name << endl;
    assert( !isReset );
    // TODO: disable the checks on release.
    activeSet.reset();
    std::fill( freeML_.begin(),freeML_.end(),true );
    isReset=true; isInit = false; isOptimumCpt = false; isLagrangeCpt = false;
    sotDEBUG(45) << "# Out }" << name << endl;
  }
  //    ,isReset(false),isInit(false),isOptimumCpt(false),isLagrangeCpt(false)



  /* Compute ML=J(initIr,:)*Y. */
  void Stage::
  computeInitialJY( const ActiveSet & initialIr )
  {
    if( isAllRow(initialIr) ) { computeInitialJY_allRows(); return; }
    if( initialIr.nbActive()==0 )
      {
	sotDEBUG(5) << "Initial IR empty." << std::endl;
	ML_.setZero(); return;
      }

    activeSet=initialIr;
    sotDEBUG(5) << "initIr = " << (MATLAB)(VectorXi)activeSet << std::endl;

    ML_.setZero();
    SubMatrix<MatrixXd,RowPermutation> Jact( J,activeSet );
    Block<MatrixXd> ML = ML_.topRows(activeSet.nbActive());
    sotDEBUG(15) << "Ja = " << (MATLAB)Jact << std::endl;
    ML = Jact;
    sotDEBUG(15) << "Ja = " << (MATLAB)ML_ << std::endl;
    Y.applyThisOnTheLeft( ML );

    for( unsigned int i=0;i<activeSet.nbActive();++i )
      {	freeML_[i]=false; }
    for( unsigned int i=0;i<nr;++i )
      {
	if( activeSet.isActive(i) )
	  e_( activeSet.where(i) ) = bounds[i].getBound( activeSet.whichBound(i) );
      }

    sotDEBUG(15) << "JY = " << (MATLAB)ML_ << std::endl;
    sotDEBUG(15) << "e = " << (MATLAB)e_ << std::endl;
  }
  /* Compute ML=J(:,:)*Y. */
  void Stage::
  computeInitialJY_allRows(void)
  {
    ActiveSet Ir(nr);
    for( unsigned int i=0;i<nr;++i )
      {
	if( bounds[i].getType() != Bound::BOUND_TWIN ) continue;
	int r = Ir.activeRow( i,Bound::BOUND_TWIN );
	assert( r==Ir.nbActive()-1 );
     }
    computeInitialJY( Ir );
  }

  /* TODO: the previous rank is not usefull anymore, since it correspond to the
   * actual rank of Y. */
  unsigned int Stage::
  computeInitialCOD( const unsigned int previousRank,
		     const ActiveSet & initialIr,
		     BaseY& Yinit)
  {
    /*
     * ML=J(initIr,:)*Y;
     * rank=Ir.size();  Ir=1:rank;
     * M=submatrix(ML,1:previousRank); L=submatrix(ML,previousRank+1:end);
     *
     * A=columnMajor(L)  % A==L'
     * qr(A);
     * RotationHouseHolder_list_t Yup( A );
     * Y=Y*Yup;
     *
     * for i=rank:-1:1
     *   if( L(i,i)!= 0 ) break;
     *     nullifyLineDeficient( i );
     */

    assert( isReset&&(!isInit ) );
    isInit=true; isReset=false;

    sotDEBUG(5) << "J = " << (MATLAB)J << std::endl;

    /* Compute ML=J(initIr,:)*Y. */
    computeInitialJY(initialIr);

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
    if( (L.cols() == 0)||(L.rows()==0) )
    {
      sotDEBUG(5) << "col size of L is null, skip the end of initialization" << std::endl;
      L.setRowIndices(VectorXi());
      L.setColIndices(VectorXi());
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

    /* Y=Y*Yup; */
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
     * Jd = L.row(r);
     * foreach i in rank:-1:1
     *   if( Jd(i)==0 ) continue;
     *     gr= GR(L(Ir(i),i),Jd(i),i,r );
     *     L=gr*L;
     *     W=W*gr';
     * Ir >> r;
     * In << r;
     */
    const Index r = (in_r<0)?row-1:in_r;
    for( Index i=r-1;i>=0;--i )
      {
	//if( std::abs(L(row,i))<EPSILON ) continue; //{ ML_(Irn(row),L.getColIndices()(i)) = 0; continue; } // PSEUDOZEROS
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

  /* --- FREEZE ------------------------------------------------------------- */
  /* --- FREEZE ------------------------------------------------------------- */
  /* --- FREEZE ------------------------------------------------------------- */


  /* Freeze the equalities constraint whose lagrange multipliers are non zero.
   * In case of slack variables (when the arg is true), also modify the
   * objective function to an reachable value.
   */
  void Stage::
  freezeSlacks( const bool & slacks )
  {
    assert(isLagrangeCpt);
    sotDEBUG(55) << "# In { " << name << endl;

    EI_FOREACH( i,lambda )
      {
	if( std::abs(lambda(i,0))>EPSILON )
	  {
	    const Index il = W.getRowIndices()(i);
	    const unsigned int cstref = activeSet.whichConstraint(il);
	    Bound::bound_t btype = activeSet.whichBound(cstref);
	    assert( (btype!=Bound::BOUND_NONE)&&(btype!=Bound::BOUND_DOUBLE) );

	    if( btype!=Bound::BOUND_TWIN )
	      {
		activeSet.freeze(cstref);
		if(!slacks)
		  { sotDEBUG(5)<<"Freeze cst "<<name<<":"<<cstref <<"."<<endl; }
	      }
	    /* Modify the bound when slack is l>0. */
	    if( slacks )
	      {
		e_(il) += lambda(i,0);
		sotDEBUG(5)<<"Freeze cst "<<name<<":"<<cstref<<" to "<<e(i,0)<< endl;
	      }
	  }
      }
    sotDEBUG(55) << "# Out } " << name << endl;
  }




  /* --- DOWNDATE ----------------------------------------------------------- */
  /* --- DOWNDATE ----------------------------------------------------------- */
  /* --- DOWNDATE ----------------------------------------------------------- */

  // Return true if the rank re-increase operated at the current stage.
  bool Stage::
  downdate( const unsigned int position,
	    GivensSequence & Ydown )
  {
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
    assert( isInit ); isLagrangeCpt=false; isOptimumCpt=false;
    //sotDEBUGPRIOR(+45);
    sotDEBUG(5) << " --- DOWNDATE ----------------------------" << std::endl;
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
    assert( isInit ); isLagrangeCpt=false; isOptimumCpt=false;
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

    /* Check if one of the M's grown. */
    //TODO: search from the end: for( Index i=sizeN()-1;i>=0;--i )
    for( Index i=0;i<sizeN();++i )
      {
	if( std::abs(ML_(Irn(i),sizeM)) > EPSILON )
	  {
	    /* Remove all the non-zero compononent of ML(i+1:end,sizeM). */
	    sotDEBUG(5) << "Found a non zero at "<<i << std::endl;
	    Block<MatrixXd> ML(ML_,0,0,nr,sizeM+1);
	    for( Index j=i+1;j<sizeN();++j )
	      {
		// if( std::abs(ML_(Irn(j),sizeM))<=EPSILON )
		//   { sotDEBUG(5) << "continue..." << endl; continue; }  // PSEUDOZERO
		Givensd G1;
		G1.makeGivens(ML_(Irn(i),sizeM),ML_(Irn(j),sizeM));
		ML.applyOnTheLeft( Irn(i),Irn(j),G1.transpose());
		W_.applyOnTheRight( Irn(i),Irn(j),G1);
		assert( std::abs(ML_(Irn(j),sizeM))<EPSILON*EPSILON ); ML_(Irn(j),sizeM)=0; // PSEUDOZERO
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
    for( int i=0;i<sizeN();++i ) { ML_(Irn(i),sizeM)=0.; } // PSEUDOZEROS
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
	/* Compare ||rest||^2 = 1-(1-x)^2 ~ 2x, with x=Wrc. */
	if( std::abs(W(row,col)-1)< EPSILON*EPSILON/2 ) break;

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
    sotDEBUG(5) << " --- UPDATE ----------------------------" << std::endl;
    assert( isInit ); isLagrangeCpt=false; isOptimumCpt=false;
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

    assert( (cst.second==Bound::BOUND_SUP)||(cst.second==Bound::BOUND_INF) );
    double sign = (cst.second==Bound::BOUND_SUP)?+1:-1;
    sotDEBUG(5) << "cst=" << cst.first << " bound=" << cst.second << "  (" << sign << ")." << endl;
    e_(wrowup) = sign*bounds[cst.first].getBound(cst.second);
    sotDEBUG(5) << "bound="<<bounds[cst.first].getBound(cst.second) <<  endl;
    RowML JupY = ML_.row(wcolup);
    JupY = sign*J.row(cst.first); Y.applyThisOnTheLeft(JupY);
    double norm2=0; int rankJ=sizeM;
    for( Index i=nc-1;i>=sizeM;--i )
      {
	norm2+=JupY(i)*JupY(i);
	if( norm2>EPSILON*EPSILON )
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
    assert( isInit ); isLagrangeCpt=false; isOptimumCpt=false;
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
	if( sizeL>0 ) M.pushColBack( L.popColFront() );
	else M.pushColBack(sizeM);
	sizeM++;
      }

    sotDEBUG(5) << "M = " << (MATLAB)M << endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << endl;

#ifdef SOTH_DEBUG
    {
      MatrixXd Yex(nc,nc); Yex.setIdentity(); Ydown.applyThisOnTheLeft(Yex);
      sotDEBUG(5) << (MATLAB)Yex << endl;
    }
#endif
    for( unsigned int i=0;i<sizeA();++i )
      {
	RowL MLi = rowML(i);
	sotDEBUG(5) << "MLib = " << (MATLAB)MLi << endl;
	Ydown.applyThisOnTheLeftReduced(MLi);
	sotDEBUG(5) << "MLia = " << (MATLAB)MLi << endl;
      }
    sotDEBUG(5) << "M = " << (MATLAB)M << endl;
    sotDEBUG(5) << "L = " << (MATLAB)L << endl;

    // sizeM already increased, so sM+sL is the last col of the Hessenberg.
    if( sizeM+sizeL<=decreasePreviousRank )
      { // L increased a column.
	if( sizeL>0 )L.pushColBack( sizeM+sizeL-1 );
      }
    else if(! defDone )
      { // rank decrease ongoing...
	const int rdef = decreasePreviousRank-sizeM;
	assert( rdef<sizeL );
	nullifyLineDeficient( rdef,rdef );
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

  /* --- DIRECT ------------------------------------------------------------- */
  /* Zu=Linv*(Ui'*ei-Mi*Yu(1:rai_1,1)); */
  void Stage::computeSolution( const VectorXd& Ytu, VectorXd& Ytdu, bool init )
  {
    assert( isInit );
    if (sizeL==0)
    {
      sotDEBUG(10) << "size of L is 0, skipping solve" << std::endl;
      return;
    }
    sotDEBUG(5) << "e = " << (MATLAB)e << std::endl;

    const VectorBlock<VectorXd> ulprec = Ytu.segment( sizeM,sizeL );
    const VectorBlock<VectorXd> umprec = Ytu.head( sizeM );
    const VectorBlock<VectorXd> dum = Ytdu.head( sizeM );
    VectorBlock<VectorXd> We = Ytdu.segment( sizeM,sizeL );
    const SubMatrixXd Wr( W_,W.getRowIndices(),Ir );
    const SubMatrixXd Mr( ML_,Ir,M.getColIndices() );

    sotDEBUG(25) << "Ytu = " << (MATLAB)Ytu << std::endl;
    sotDEBUG(25) << "Ytdu = " << (MATLAB)Ytdu << std::endl;
    sotDEBUG(45) << "Wr = " << (MATLAB)Wr << std::endl;

    /* TODO: when L0 is full rank, a permuation of e should be enough (W=Id). */
    We = Wr.transpose()*e;
    sotDEBUG(25) << "Wre = " << (MATLAB)We << std::endl;
    if(! init )
      {
	We.noalias() -= getLtri()*ulprec;
      }
    sotDEBUG(25) << "Wre_Lu = " << (MATLAB)We << std::endl;
    if( sizeM >0 )
      {
	We.noalias() -= Mr*(umprec+dum); /* TODO: this sum u+du could be done
					  * only once, while it is done at each
					  *  stage now. */
      }
    sotDEBUG(5) << "Wre_Lu_Mru_Mrdu = " << (MATLAB)We << std::endl;

    soth::solveInPlaceWithLowerTriangular(L,We);
    sotDEBUG(5) << "LiWrde = " << (MATLAB)We << std::endl;
    sotDEBUG(45) << "Ytdu = " << (MATLAB)Ytdu << std::endl;
  }

  /* --- INDIRECT ----------------------------------------------------------- */

  /* err = Ju-e = W [M L 0] Y^u - e
   * where MLYtu has already been computed.
   */
  void Stage::
  computeErrorFromJu(const VectorXd& MLYtu, VectorXd& err) const
  {
    assert(MLYtu.size() == sizeA());
    err.noalias() =  W*MLYtu-e; // DEBUG  e-W*MLYtu
  }

  /* Compute W'Ju = MLYtu = M*Ytu.head + L*Ytu.tail.
   */
  void Stage::
  computeMLYtu( const VectorXd& Ytu,VectorXd& MLYtu ) const
  {
    assert(Ytu.size() == nc);
    MLYtu.noalias() = M*Ytu.head(sizeM);
    MLYtu.tail(sizeL).noalias()
      += getLtri()*Ytu.segment(sizeM, sizeL);

    sotDEBUG(45) << "Ytu = " << (MATLAB)Ytu << endl;
    sotDEBUG(45) << "M = " << (MATLAB)M << endl;
    sotDEBUG(45) << "L = " << (MATLAB)(MatrixXd)getLtri() << endl;
    sotDEBUG(45) << "MLYtu = " << (MATLAB)MLYtu << endl;

  }

  /* err = Ju-e = W [M L 0] Y^u - e
   */
  void Stage::
  computeError(const VectorXd& Ytu, VectorXd& err) const
  {
    assert(Ytu.size() == nc);
    assert( isInit );

    // TODO: temporary allocation in this expression. Manage it?
    VectorXd MLYtu; computeMLYtu(Ytu,MLYtu);
    computeErrorFromJu(MLYtu,err);
  }

  /* Compute the error from scratch, and stored it in lambda. */
  void Stage::
  computeError(const VectorXd& Ytu)
  {
    assert( isInit );
    VectorXd tmp;
    computeError(Ytu,tmp);

    const Indirect & Ie = e.getRowIndices();
    lambda.setRowIndices( Ie );
    TRANSFERT_IN_SUBVECTOR(tmp,lambda);
    isLagrangeCpt =true;
  }

  /* Compute the error from already compute MLYtu, and stored it in lambda. */
  void Stage::
  computeErrorFromJu(const VectorXd& MLYtu)
  {
    VectorXd tmp;
    computeErrorFromJu(MLYtu,tmp);

    const Indirect & Ie = e.getRowIndices();
    lambda.setRowIndices( Ie );
    TRANSFERT_IN_SUBVECTOR( tmp,lambda );
    isLagrangeCpt =true;
  }



  /* Compute J' (Ju-e) in Y base:
   * Ytrho = -Y'J'(Ju-e) = [ M L ]' ( W'e - [M L] Ytu ). */
  void Stage::
  computeRho(const VectorXd& Ytu, VectorXd& Ytrho, bool inLambda )
  {
    assert( isInit );
    assert(Ytu.size() == nc);

    //TODO : manage temporary memory ?
    VectorXd MLYtu; computeMLYtu(Ytu,MLYtu);
    if( inLambda )
      {
	/* Compute lambda = e-W*MLYtu by the way. */
	computeErrorFromJu(MLYtu);
      }
    sotDEBUG(5) << "WtJu = " << (MATLAB)MLYtu << endl;

    /* MLYtu := W'e - [ML] Yt u . */
    MLYtu *= -1; MLYtu.noalias() += W.transpose()*e;
    sotDEBUG(5) << "Wte_Ju = " << (MATLAB)MLYtu << endl;

    /* TODO: W'(e-Ju) is null on the Ir part.
     * ... maybe not for u=u0+du, with non null u0. */

    /* Ytrho := [M L]' * MLYtu = [ M L ]' ( W'e - [M L] Ytu ). */
    Ytrho.resize( nc ); /* TODO: should be resize only once at the construction of the HCOD. */
    Ytrho.head(sizeM).noalias() = M.transpose()*MLYtu;
    Ytrho.segment(sizeM,sizeL).noalias()
      = getLtri().transpose()*MLYtu.tail(sizeL);
    Ytrho.tail(nc-sizeM-sizeL).setZero(); /* TODO: nobody ever will access the tail, could be neglected. */
  }


   /** input: rho_under_i = {ro_1, ..., ro_i}
    * on return:
    * lambda_i =  Wr_i*L_i^{-T}*rho_i
    * rho_under_{i-1} = rho_under_{i-1} + Mr_i^T*L_i^{-T}*rho_i
//???    * rho_i = L_i^{-T}*rho_i (should not be useful).
    */
  void Stage::
  computeLagrangeMultipliers( VectorXd& rho, VectorXd& l ) const
  {

    if( 1 )
      {
	assert( isInit );
	assert( rho.rows() == nc );

	if( sizeL==0 )
	  {	l.resize(sizeA()); l.setZero(); return; }
	VectorBlock<VectorXd> rho_i = rho.segment(sizeM,sizeL);
	sotDEBUG(5) << "rho = " << (MATLAB)rho_i << endl;

	solveInPlaceWithUpperTriangular(L.transpose(), rho_i);
	sotDEBUG(5) << "Lirho = " << (MATLAB)rho_i << endl;

	l.noalias() = W.rightCols(sizeL) * rho_i;
	if( sizeM>0 )
	  {
	    VectorBlock<VectorXd> rho_under = rho.head(sizeM);
	    rho_under.noalias() -= M.bottomRows(sizeL).transpose()*rho_i;
	  }
      }
    else // DEBUG!!
      {
	MatrixXd J_(nr,nc); VectorXd e_(nr);
	SubMatrix<MatrixXd,RowPermutation> Ja = Jactive(J_);
	SubVectorXd ea = eactive(e_);
	l = Ja*rho-ea;
      }
  }

  void Stage::
  computeLagrangeMultipliers( VectorXd& rho )
  {
    assert( isInit );
    // TODO: deal with the temporary.
    VectorXd ltmp; computeLagrangeMultipliers(rho,ltmp);
    lambda.setRowIndices( e.getRowIndices() );
    TRANSFERT_IN_SUBVECTOR(ltmp,lambda);
    sotDEBUG(1) << "l = " << (MATLAB)lambda << endl;
    isLagrangeCpt=true;
  }

  void Stage::
  transfertInSubVector( const VectorXd& tmp, VectorXd& rec,const Indirect& idx )
  {
    EI_FOREACH( i,idx )
      {
	rec(idx(i)) = tmp(i);
      }
  }


  /* --- BOUND -------------------------------------------------------------- */

  /* Return true if the previous tau is correct (ie maxlocal(tau)>taumax */
  bool Stage::checkBound( const VectorXd& u,const VectorXd& du,
			  ConstraintRef& cstmax, double& taumax )
  {
    bool res = true;
    for( unsigned int i=0;i<nr;++i )
      {
	if( activeSet.isActive(i) ) continue;
	assert( bounds[i].getType()!=Bound::BOUND_TWIN );

	/* This has already been computed and could be avoided... TODO. */
	double val = J.row(i)*u;
	double dval = J.row(i)*du;
	const Bound & b = bounds[i];
	sotDEBUG(5) << "bound = " << b << endl;
	sotDEBUG(5) << "Ju = " << val << "  --  Jdu = " << dval  << " -- Jupdu = " << val+dval << endl;


	Bound::bound_t btype = b.check(val,EPSILON);
	if( btype!=Bound::BOUND_NONE )
	  {
	    assert( (btype==Bound::BOUND_INF)||(btype==Bound::BOUND_SUP) );
	    sotDEBUG(5) << "Violation Ju at " <<name <<" " << ((btype==Bound::BOUND_INF)?"-":"+")<<i << std::endl;
	    taumax=0; cstmax = std::make_pair(i,btype);
	    return false;
	  }
	else
	  {
	    btype = b.check(val+dval,EPSILON);//DEBUG
	    if( btype!=Bound::BOUND_NONE )
	      {
		assert( (btype==Bound::BOUND_INF)||(btype==Bound::BOUND_SUP) );
		sotDEBUG(5) << "Violation at " <<name <<" "<< ((btype==Bound::BOUND_INF)?"-":"+")<<i << std::endl;

		const double & bval = b.getBound(btype);
		double btau = (bval-val)/dval;
		assert(btau>=0); assert(btau<1);
		if( btau<taumax )
		  {
		    sotDEBUG(1) << "Max violation (tau="<<btau<<") at "<<name <<" "
				<< ((btype==Bound::BOUND_INF)?"-":"+")<<i << std::endl;
		    res=false;
		    taumax=btau; cstmax = std::make_pair(i,btype);
		  }
	      }
	  }
      }
    return res;
  }

  bool Stage::checkBound( const VectorXd& u,const VectorXd& du,
			  ConstraintRef* cstptr, double* tauptr )
  {
    if( (tauptr==NULL)&&(cstptr==NULL) )
      { double tau=1; ConstraintRef cst; return checkBound(u,du,cst,tau); }
    else if( (tauptr!=NULL)&&(cstptr==NULL) )
      { ConstraintRef cst; return checkBound(u,du,cst,*tauptr); }
    else if( (tauptr==NULL)&&(cstptr!=NULL) )
      { double tau=1; return checkBound(u,du,*cstptr,tau); }
    else if( (tauptr!=NULL)&&(cstptr!=NULL) )
      { return checkBound(u,du,*cstptr,*tauptr); }


  }

  bool Stage:: // TODO: Ytu could be passed instead of u.
  maxLambda( const VectorXd& u, double & lmax,unsigned int& row ) const
  {
    /* TODO: unactive the search for TWINS. */

    bool res=false;
    EI_FOREACH( i,lambda )
      {
	const unsigned int cstref = activeSet.whichConstraint(W.getRowIndices()(i));
	Bound::bound_t btype = activeSet.whichBound(cstref);
	assert( (btype!=Bound::BOUND_NONE)&&(btype!=Bound::BOUND_DOUBLE) );
	switch( btype )
	  {
	  case Bound::BOUND_TWIN:
	    break; // Nothing to do.
	  case Bound::BOUND_SUP:
	    sotDEBUG(5) << name<<": row"<<i<<", cst"<<which(i) << ": l=" << lambda(i,0) << endl;
	    if( -lambda(i,0)>lmax )
	      {
		double Ju = J.row(cstref)*u;
		if( Ju<=bounds[cstref].getBound( Bound::BOUND_SUP )+EPSILON )
		  {
		    res=true;
		    lmax=-lambda(i,0);
		    row=i;
		  }
	      }
	    break;
	  case Bound::BOUND_INF:
	    sotDEBUG(5) << name<<": row"<<i<<", cst"<<which(i) << ": l=" << lambda(i,0) << endl;
	    if( -lambda(i,0)>lmax ) // TODO: change the sign of the bound-inf cst.
	      {
		double Ju = J.row(cstref)*u;
		if( bounds[cstref].getBound( Bound::BOUND_INF )-EPSILON<=Ju )
		  {
		    res=true;
		    lmax=-lambda(i,0);
		    row=i;
		  }
	      }
	    break;
	  }
      }
    return res;
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
  {
    if( r<sizeN() ) return sizeM;
    else return std::min( sizeM+r-sizeN()+1,nc );
    //return (r<sizeN())?sizeM:sizeM+r-sizeN()+1;
 }

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
    bool res; double norm=0;
    if( sizeA() ) { norm=(Jrec-Ja).norm(); res = (norm<=10*EPSILON); }
    else res = ( (Jrec.cols()==Ja.cols())&&(Jrec.rows()==Jrec.rows()) );
    sotDEBUG(5) <<"% J: Recomposition  " << ((res)?"OK":"wrong")
		<< " (||.||="<<norm<<")."<< std::endl;

    VectorXd ea_;
    bool vres=true; // DEBUG: the test is wrong when the stage has been freezed.
    // if( sizeA() ) vres = (e-eactive(ea_)).norm()<=EPSILON;
    // else vres = (e.size()==eactive(ea_).size() );

    sotDEBUG(15) << "e="<<(MATLAB)e << endl;
    sotDEBUG(15) << "ea="<<(MATLAB)eactive(ea_) << endl;
    sotDEBUG(5) <<"% e: Recomposition  " << ((vres)?"OK.":"wrong.") << std::endl;

    return res&&vres;
  }

  Stage::Index Stage::
  where( unsigned int cst ) const
  {
    Index ref = activeSet.where(cst);
    const Indirect & Idx = W.getRowIndices();
    for( Index i=0;i<nr;++i )
      { if( Idx(i)==ref )  return i; }
  }
  Stage::ConstraintRef Stage::
  which( unsigned int row ) const
  {
    assert( row<sizeA() );
    ConstraintRef res;
    res.first = activeSet.whichConstraint( W.getRowIndices()(row) );
    res.second = activeSet.whichBound( res.first );
    return res;
  }
  bool Stage::
  isActive( unsigned int cst ) const
  {
    assert( cst<nr );
    return activeSet.isActive(cst);
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
	    J_.row(activeSet.where(i)) = activeSet.sign(i)*J.row(i);
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
	    e_(activeSet.where(i)) = activeSet.sign(i)*bounds[i].getBound(activeSet.whichBound(i));
	  }
      }

    return  SubVectorXd(e_,W.getRowIndices());
  }

  void Stage::
  show( std::ostream& os, unsigned int stageRef, bool check ) const
  {
    sotDEBUGIN(5);

    if( sotDEBUG_ENABLE(55) && sotDEBUGFLOW.outputbuffer.good() ) activeSet.disp( sotDEBUGFLOW.outputbuffer );
    os << "sa{"<<stageRef<<"} = " << (MATLAB)sizeA() << endl;
    os << "r{"<<stageRef<<"} = " << (MATLAB)sizeL << endl;
    os << "sn{"<<stageRef<<"} = " << (MATLAB)sizeN() << endl;
    os << "sm{"<<stageRef<<"} = " << (MATLAB)sizeM << endl;

    MatrixXd J_(nr,nc); J_.setConstant(-1.11111);
    VectorXd e_(nr); e_.setConstant(-1.11111);
    for( unsigned int i=0;i<nr;++i )
      {
	if( activeSet.isActive(i) )
	  {
	    double sign = (activeSet.whichBound(i)==Bound::BOUND_INF)?-1:+1;
	    J_.row(activeSet.where(i)) = sign*J.row(i);
	    //DEBUGe_(activeSet.where(i)) = sign*bounds[i].getBound(activeSet.whichBound(i));
	    sotDEBUG(55) << "where(" << i << ") = " << activeSet.where(i) << endl;
	  }
      }

    SubMatrix<MatrixXd,RowPermutation> Ja(J_,W.getRowIndices());
    SubVectorXd ea(e_,W.getRowIndices());

    sotDEBUG(55) << "Iw1{"<<stageRef<<"} = " << (MATLAB)W.getRowIndices() << std::endl;
    sotDEBUG(5) << "Iw2{"<<stageRef<<"} = " << (MATLAB)W.getColIndices() << std::endl;
    sotDEBUG(55) << "Ie{"<<stageRef<<"} = " << (MATLAB)e.getRowIndices() << std::endl;
    sotDEBUG(55) << "Il{"<<stageRef<<"} = " << (MATLAB)L.getRowIndices() << std::endl;
    sotDEBUG(55) << "J_{"<<stageRef<<"} = " << (MATLAB)J_ << std::endl;
    sotDEBUG(25) << "ML_{"<<stageRef<<"} = " << (MATLAB)ML_ << std::endl;
    sotDEBUG(5) << "erec{"<<stageRef<<"} = " << (MATLAB)ea << std::endl;

    os << "a{"<<stageRef<<"} = " << (MATLAB)(Indirect)activeSet << std::endl;
    os << "J{"<<stageRef<<"} = " << (MATLAB)Ja << std::endl;
    os << "e{"<<stageRef<<"} = " << (MATLAB)e << std::endl;
    os << "W{"<<stageRef<<"} = " << (MATLAB)W << std::endl;
    os << "M{"<<stageRef<<"} = " << (MATLAB)M << std::endl;
    os << "L{"<<stageRef<<"} = " << (MATLAB)L << std::endl;

    if( check )
    {
      MatrixXd Jrec; recompose(Jrec);
      if (Jrec.rows()>0)
      {
        sotDEBUG(5) << "Jrec="<<(MATLAB)Jrec << endl;
        if((Jrec-Ja).norm()>1e-6) os << "Jrec{"<<stageRef<<"} = " << (MATLAB)Jrec << std::endl;
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
    if( isLagrangeCpt )
      {
	os << "lag{"<<stageRef<<"} = " << (MATLAB)lambda << std::endl;
      }
  }

  void Stage::
  showActiveSet( std::ostream& os ) const
  {
    os << " [ ";
    for( unsigned int r=0;r<nr;++r )
      {
	if(! activeSet.isActive(r) ) continue;
	if( bounds[r].getType() ==  Bound::BOUND_DOUBLE )
	  if( activeSet.whichBound(r) == Bound::BOUND_INF ) os << "-";
	  else if( activeSet.whichBound(r) == Bound::BOUND_SUP ) os << "+";
	os <<r << " ";
      }
    os << " ]";
  }


  std::ostream& operator<<( std::ostream&os,const Stage::ConstraintRef& cst )
  {
    switch( cst.second )
      {
      case Bound::BOUND_INF: os << "-"; break;
      case Bound::BOUND_SUP: os << "+"; break;
      case Bound::BOUND_DOUBLE: os << "+/-"; break;
      case Bound::BOUND_TWIN: os << "="; break;
      case Bound::BOUND_NONE: os << "(o)"; break;
      }
    return os << cst.first;
  }

  ActiveSet Stage::_allRows(0);
  double Stage::EPSILON = 1e-6;

}; // namespace soth
