#ifndef __SOTH_STAGE__
#define __SOTH_STAGE__

#include <Eigen/Core>
#include <Eigen/Jacobi>
#include <list>
namespace Eigen
{
  #include "soth/SubMatrix.h"
}
#include "soth/Algebra.h"
#include "soth/BaseY.hpp"
#include "soth/Bound.hpp"
#include "soth/ActiveSet.hpp"

namespace soth
{


  /* --- STAGE -------------------------------------------------------------- */
  /* --- STAGE -------------------------------------------------------------- */
  /* --- STAGE -------------------------------------------------------------- */

  class Stage
  {
  public:
    typedef MatrixXd::Index Index;
    typedef SubMatrix<MatrixXd> SubMatrixXd;
    typedef SubMatrix<VectorXd,RowPermutation> SubVectorXd;
    typedef SubMatrixXd::RowIndices Indirect;

    typedef VectorBlock<MatrixXd::RowXpr> RowL;
    typedef MatrixXd::RowXpr RowML;

    typedef std::pair<Index,Bound::bound_t> ConstraintRef;

    typedef PlanarRotation<double> Givensd;
    typedef std::list<Givensd> givensd_sequence_t;

  protected:

    const MatrixXd & J;
    const bound_vector_t & bounds;
    BaseY & Y;

    unsigned int nr,nc; // nr=nbCols(J), nc=nbRows(J).
    ActiveSet activeSet;

    MatrixXd W_;
    MatrixXd ML_;
    VectorXd e_;

    SubMatrixXd M,L;
    SubVectorXd e;

    bool isWIdenty;
    SubMatrixXd W;

    //SubMatrixXd L0sq;
    //TriSubMatrixXd L0; // L0 = tri(L0sq)

    //TriMatrixXd Ldamp;

    // fullRankRows = Ir. defRankRows = In.
    const Indirect& Ir, &Irn; // Ir = L0sq.indirect1() -- Irn = M.indirect1()

    unsigned int sizeM,sizeL; // sizeL = card(Ir).
    unsigned int sizeN; // sizeN = card(In) = sizeA-sizeL.

    /* W = W_( :,[In Ir] ).
     * M = ML_( [In Ir],0:sizeM-1 ).
     * L = ML_( Ir,sizeM:sizeM+sizeL-1  ).
     * Lo = [0;L] = ML_( [In Ir],sizeM:sizeM+sizeL-1  ).
     *
     * W_*ML_ = W*[M [ zeros(In.size(),rank); L ] ]
     */

  public:

    Stage( const MatrixXd & J, const bound_vector_t & bounds, BaseY& Y  );

    /* --- INIT ------------------------------------------------------------- */

    /* Return the rank of the current COD = previousRank+size(L). */
    unsigned int computeInitialCOD( const unsigned int previousRank,
				    const ActiveSet & initialIr );
    /*
      ML=J(initIr,:)*Y;
      rank=Ir.size();  Ir=1:rank;
      M=submatrix(ML,1:previousRank); L=submatrix(ML,previousRank+1:end);

      A=columnMajor(L)  // A==L'
      qr(A);
      RotationHouseHolder_list_t Yup( A );
      Y=Y*Yup;

      for i=rank:-1:1
      if( L(i,i)!= 0 ) break;
      nullifyLineDeficient( i );

    */

  protected:
    void nullifyLineDeficient( const Index row, const Index in_r );
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

    void computeInitalJY( const ActiveSet & initialIr );
    void computeInitalJY_allRows(void);

    /* --- DOWN ------------------------------------------------------------- */
  public:
    // Return true if the rank re-increase operated at the current stage.
    bool downdate( const unsigned int position,
		   givensd_sequence_t & Ydown );
    /*
      gr = Neutralize row <position> in W <position>
      L=gr'*L
      bool res;
      if( L(In.last(),0) == 0
      // Rank deficience: regularize Hessenberg
      Ydown = regularizeHessenberg
      res=false;
      else
      // No rank dec: quit
      Ir << In.pop();
      res=true;
      Ir >> r; In >> r;
      Unused << r;
      return res;
    */


    // Return true if the rank decrease operated at the current stage.
    bool propagateDowndate( givensd_sequence_t & Ydown,
			    bool decreasePreviousRank );
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

    void regularizeHessenberg( givensd_sequence_t & Ydown );
    /*
      for i=i0:rank-1
      gr = GR( L(Ir(i),i),L(Ir(i),i+1),i,i+1 );
      L = L*GR;
      Ydown.push_back( gr );
    */
  protected:
    void removeInW( const  unsigned int position );

    /* --- UPD -------------------------------------------------------------- */
  public:
    /* TODO
       if stage(sup).update( Jup,Yup )
       while stage(sup++).propagateUpdate( Yup,false
       else ..

    */

    // Return true if the rank re-decrease operated at the current stage.
    bool update( const VectorXd & Jup,
		 givensd_sequence_t & Yup );
    /*
      Inew = Unused.pop();
      Row JupY = row(Inew);
      JupU = Jup*Y;
      double norm2=0; double rankJ=0;
      for i=n:-1:1
      norm2+=JupY(i)^2;
      if norm2!=0 rankJ=i; break;

      Ir << Inew
      W(Inew,Inew)=1;
      if rankJ>sizeM+rank
      // Rank increase
      for i=rankJ:-1:sizeM+rank+1
      Gr gr; gr.init( JupY,i,i-1,0 ); prod(JupY,gr);
      Yup.push_back( gr );
      return false;
      else
      // No rank increase;
      nullifyLineDeficient(Inew);
      return true;
    */


    // Return true if the rank decrease operated at the current stage.
    //bool propagateUpdate( RotationGiven_list_t & Yup,
    //bool decreasePreviousRank );
    /*
      ML=ML*Yup;
      increaseSizeM();
      if( decreasePreviousRank ) return true
      for i=1:rank
      if L(i,i) == 0
      // Rank re-decrease
      nullifyLineDeficient( i );
      regularizeHessenberg( Yup,i+1 );
      return true;
      return false;
    */


    /* --- SOLVE ------------------------------------------------------------ */

    /* Solve in the Y space. The solution has then to be multiply by Y: u = Y*Yu. */
    void solve( VectorXd& Yu );
    //void solveTranspose( const VectorXd & e, VectorXd res );

    //void damp( const double & damping );


    /* --- CHECK ------------------------------------------------------------ */
    /* WMLY = [ W*M W(:,1:rank)*L zeros(sizeA,nc-sizeM-sizeL) ]*Y' */
    void recompose( MatrixXd& WMLY );
    void show( std::ostream& os, unsigned int stageRef, bool check=false );

  public:
    /* --- ACCESSORS --- */

    // SubMatrixXd M();
    // SubMatrixXd L();
    // SubMatrixXd Lo();

    // const_SubMatrixXd M() const ;
    // const_SubMatrixXd L() const ;
    // const_SubMatrixXd Lo() const ;

    RowL rowL0( const unsigned int r );
    RowML rowMrL0( const unsigned int r );
    RowL rowML( const unsigned int r );
    unsigned int rowSize( const unsigned int r );


    int sizeA( void ) const { return activeSet.nbActive(); }

    // SubRowXd row( const unsigned int r );
    // SubRowXd rown( const unsigned int r ); // row rank def
    // SubRowXd rowf( const unsigned int r ); // row rank full

    /* --- MODIFIORS --- */
    // void increaseSizeM();
    // void decreaseSizeM();

  public:
    static ActiveSet& allRows() { return _allRows; }
    static double EPSILON;
  protected:
    static ActiveSet _allRows;
    bool isAllRow( const ActiveSet& idx ) { return (&idx == &_allRows); }

  };

}; // namespace soth


#endif // #ifndef __SOTH_STAGE__
