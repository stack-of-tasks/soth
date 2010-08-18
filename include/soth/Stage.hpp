#ifndef __SOTH_STAGE__
#define __SOTH_STAGE__

#include <Eigen/Core>
#include <list>
#include <string>
#include "soth/SubMatrix.hpp"
#include "soth/solvers.hpp"
#include "soth/Algebra.hpp"
#include "soth/BaseY.hpp"
#include "soth/Bound.hpp"
#include "soth/ActiveSet.hpp"
#include "soth/Givens.hpp"

namespace soth
{


  /* --- STAGE -------------------------------------------------------------- */
  /* --- STAGE -------------------------------------------------------------- */
  /* --- STAGE -------------------------------------------------------------- */

  class Stage
  {
  public:
    typedef MatrixXd::Index Index;
    typedef SubMatrix<MatrixXd>::RowIndices Indirect;

    typedef SubMatrix<MatrixXd> SubMatrixXd; // TODO: Mv this in SubMatrix
    typedef SubMatrix<VectorXd,RowPermutation> SubVectorXd; // TODO: Mv this in SubMatrix
    typedef SubMatrixXd const_SubMatrixXd; // TODO: Mv this in SubMatrix
    typedef SubVectorXd const_SubVectorXd; // TODO: this too
    typedef TriangularView<SubMatrixXd,Lower> TriSubMatrixXd; // TODO: usefull
    typedef TriangularView<const_SubMatrixXd,Lower> const_TriSubMatrixXd;

    typedef VectorBlock<MatrixXd::RowXpr> RowL;
    typedef MatrixXd::RowXpr RowML;

    typedef PlanarRotation<double> Givensd; // TODO: remove this type

  protected:

    const MatrixXd & J;
    const bound_vector_t & bounds;
    const BaseY & Y;

    unsigned int nr,nc; // nr=nbCols(J), nc=nbRows(J).

    //ActiveSet activeSet;
    std::vector<bool> freeML_;

    MatrixXd W_;
    MatrixXd ML_;
    VectorXd e_;
    VectorXd lambda_;

    SubMatrixXd M,L;

    bool isWIdenty;
    SubMatrixXd W;

    //TriMatrixXd Ldamp;

    /* fullRankRows = Ir. defRankRows = In.
     * Ir = L.indirect1() -- Irn = M.indirect1(). */
    Indirect& Ir, &Irn, &Iw, &Im, &Il;
    /* sizeL = card(Ir). sizeM = previousRank. */
    unsigned int sizeM,sizeL;


    SubMatrixXd Wr,Mr;
    SubVectorXd e,lambda;


    SubActiveSet<ActiveSet,Indirect> activeSet;

    /* W = W_( :,[In Ir] ).
     * M = ML_( [In Ir],0:sizeM-1 ).
     * L = ML_( Ir,sizeM:sizeM+sizeL-1  ).
     * Lo = [0;L] = ML_( [In Ir],sizeM:sizeM+sizeL-1  ).
     *
     * W_*ML_ = W*[M [ zeros(In.size(),rank); L ] ]
     */

    /* Check if the stage has been reset, initialized, if the optimum
     * has been computed, and if the lagrange multipliers have been
     * computed. */
    bool isReset,isInit,isOptimumCpt,isLagrangeCpt;
  public:

    Stage( const MatrixXd & J, const bound_vector_t & bounds, BaseY& Y  );

    /* --- INIT ------------------------------------------------------------- */

    void reset();
    /* Return the rank of the current COD = previousRank+size(L).
     * Give a non-const ref on Y so that it is possible to modify it.
     */
    unsigned int computeInitialCOD( const unsigned int previousRank,
				    const ActiveSet & initialIr,
				    BaseY & Yinit );

  protected:
    void nullifyLineDeficient( const Index row, const Index in_r );
    void computeInitialJY( const ActiveSet & initialIr );
    void computeInitialJY_allRows(void);

    /* --- DOWN ------------------------------------------------------------- */
  public:
    /* Return true if the rank re-increase operated at the current stage. */
    bool downdate( const unsigned int position,
		   GivensSequence & Ydown );
    /* Return true if the rank decrease operated at the current stage. */
    bool propagateDowndate( GivensSequence & Ydown,
			    bool decreasePreviousRank );

  protected:
    void regularizeHessenberg( GivensSequence & Ydown );
    unsigned int nullifyACrossFromW( const  unsigned int position );
    void removeARowFromL( unsigned int row );
    void removeACrossFromW( const unsigned int & row, const unsigned int & col );

    /* --- UPD -------------------------------------------------------------- */
  public:
    /* Return the rank of the line where the rank re-decrease will occurs. */
    unsigned int update( const ConstraintRef & cst, GivensSequence & Yup );
    void propagateUpdate( GivensSequence & Ydown,
			  unsigned int decreasePreviousRank );
  protected:
    void addARow( const Index & mlrowup,bool deficient=false );

    /* --- SOLVE ------------------------------------------------------------ */
  public:
    /* Solve in the Y space. The solution has then to be multiply by Y: u = Y*Yu. */
    void computeSolution( const VectorXd& Ytu, VectorXd & Ytdu, bool init );
    //void damp( const double & damping );

    /* --- MULTIPLICATORS --------------------------------------------------- */
  public:
    /* The const functions simultaneously set up the lambda member. */
    template <typename D>
    void computeError(const VectorXd& Ytu, MatrixBase<D>& err ) const;
    void computeError(const VectorXd& Ytu );
    /* computeRho(.,.,false) is const. What trick can we use to explicit that? TODO. */
    void computeRho(const VectorXd& Ytu, VectorXd& Ytrho, bool inLambda = false );
    template <typename D>
    void computeLagrangeMultipliers( VectorXd& rho, MatrixBase<D>& l ) const;
    void computeLagrangeMultipliers( VectorXd& rho );

  protected:
    void computeErrorFromJu(const VectorXd& MLYtu);
    template <typename D>
    void computeErrorFromJu(const VectorXd& Ytu, MatrixBase<D>& err) const;
    void computeMLYtu( const VectorXd& Ytu,VectorXd& MLYtu ) const;

    /* --- ACTIVE SEARCH ---------------------------------------------------- */
  public:
    /* Return true if all bounds are checked with the specified tau.  If tau is
     * specified, the step is computed btw (with tau_out <= tau_in) and the
     * constraint to update is returned.
     */
    bool checkBound( const VectorXd& u,const VectorXd& du,
		     ConstraintRef*, double* tau );
    bool checkBound( const VectorXd& u,const VectorXd& du,
		     ConstraintRef& cstmax, double& taumax );
    bool maxLambda( const VectorXd& u, double & lmax, unsigned int & row ) const;
    void freezeSlacks(const bool & slacks = true);



    /* --- CHECK ------------------------------------------------------------ */
  public:
    /* WMLY = [ W*M W(:,1:rank)*L zeros(sizeA,nc-sizeM-sizeL) ]*Y' */
    void recompose( MatrixXd& WMLY ) const;
    void show( std::ostream& os, unsigned int stageRef, bool check=false ) const;
    void showActiveSet( std::ostream& os ) const;

    /* Return a sub matrix containing the active rows of J, in the
     * same order as given by W. J_ is a matrix where th full
     * J is stored (workspace). */
    SubMatrix<MatrixXd,RowPermutation>   Jactive( MatrixXd& J_ ) const ;
    MatrixXd  Jactive() const ;
    /* Return a sub vector containing the active rows of e, in the
     * same order as given by W. */
    SubMatrix<VectorXd,RowPermutation>  eactive( VectorXd& e_ ) const;
    VectorXd  eactive() const ;

    bool testRecomposition( void ) const;
    bool testSolution( const VectorXd & solution ) const;

    /* For debug purpose, give the line of an active constraint (assert the activity). */
    Index where( unsigned int cst ) const;
    ConstraintRef which( unsigned int row ) const;
    bool isActive( unsigned int cst ) const;

  public:
    /* --- ACCESSORS --- */

    SubMatrixXd getM() { return M; }
    const_SubMatrixXd getM() const { return M; }
    SubMatrixXd getL() { return L; }
    const_SubMatrixXd getL() const { return L; }
    TriSubMatrixXd getLtri() { return L.triangularView<Lower>(); }
    const_TriSubMatrixXd getLtri() const { return L.triangularView<Lower>(); }
    SubVectorXd gete() { return e; }
    const_SubVectorXd gete() const { return e; }
    SubVectorXd getLagrangeMultipliers() { return lambda; }
    const_SubVectorXd getLagrangeMultipliers() const { return lambda; }
    template< typename D>
    void setLagrangeMultipliers(const MatrixBase<D>& l) // DEBUG
    {
      const Indirect & Ie = e.getRowIndices();
      lambda.setRowIndices( Ie );
      TRANSFERT_IN_SUBVECTOR(l,lambda);
      isLagrangeCpt =true;
    } // DEBUG

    RowL rowL0( const Index r );
    RowML rowMrL0( const Index r );
    RowL rowML( const Index r );
    unsigned int rowSize( const Index r );

    /* TODO: sizeL and sizeM should be automatically determined from the corresponding indexes. */
    unsigned int nbConstraints( void ) const { return nr; }
    unsigned int sizeA( void ) const { return activeSet.nbActive(); }
    // sizeN = card(In) = sizeA-sizeL.
    int sizeN( void ) const { assert(sizeA()-sizeL>=0);return sizeA()-sizeL; }

    Index rank() const {return sizeL;}

  public:
    static ActiveSet& allRows() { return _allRows; }
    static double EPSILON;
  protected:
    static ActiveSet _allRows;
    bool isAllRow( const ActiveSet& idx ) { return (&idx == &_allRows); }
  public: /* For debug purpose, could be remove on RELEASE. */
    std::string name;

  };

}; // namespace soth


#endif // #ifndef __SOTH_STAGE__
