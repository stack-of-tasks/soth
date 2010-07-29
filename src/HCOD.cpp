#define SOTH_DEBUG
#define SOTH_DEBUG_MODE 15
#include "soth/HCOD.hpp"
#include "soth/debug.h"

namespace soth
{

  HCOD::
  HCOD( unsigned int inSizeProblem, unsigned int nbStage )
  :
    sizeProblem(inSizeProblem)
    ,Y(sizeProblem)
    ,stages(0),initialActiveSets(0)
    ,solution(sizeProblem)
    ,du(sizeProblem),Ytu(sizeProblem),Ytdu(sizeProblem),rho(sizeProblem)
    ,isReset(false),isInit(false),isSolutionCpt(false)
  {
    stages.reserve(nbStage);
    initialActiveSets.reserve(nbStage);
  }

  /* --- SETTERS/GETTERS ---------------------------------------------------- */
  /* --- SETTERS/GETTERS ---------------------------------------------------- */
  /* --- SETTERS/GETTERS ---------------------------------------------------- */

  void HCOD::
  pushBackStage( const MatrixXd & J, const bound_vector_t & bounds )
  {
    unsigned int s = stages.size();
    stages.resize( s+1 );
    stages[s] = stage_ptr_t(new soth::Stage( J,bounds,Y ));

    initialActiveSets.resize( s+1 );
    initialActiveSets[s].resize(0);
    isInit=false;
  }
  void HCOD::
  pushBackStage( const MatrixXd & J, const bound_vector_t & bounds,const VectorXi& Ir0 )
  {
    pushBackStage(J,bounds);
    setInitialActiveSet(Ir0,stages.size()-1);
  }
  void HCOD::
  pushBackStages( const std::vector<MatrixXd> & J,
		  const std::vector<bound_vector_t> & bounds )
  {
    assert( J.size() == bounds.size() );
    for( unsigned int i=0;i<J.size();++i )
      {
	pushBackStage( J[i],bounds[i] );
      }
  }

  Stage& HCOD::
  stage( unsigned int i )
  {
    assert( i<stages.size() );
    return *stages[i];
  }
  const Stage& HCOD::
  stage( unsigned int i ) const
  {
    assert( i<stages.size() );
    return *stages[i];
  }

  void HCOD::
  setInitialActiveSet( const VectorXi& Ir0,unsigned int i )
  {
    assert( i<initialActiveSets.size() );
    initialActiveSets[i] = Ir0;
  }
  const VectorXi& HCOD::
  getInitialActiveSet( unsigned int i )
  {
    assert( i<initialActiveSets.size() );
    return initialActiveSets[i];
  }

  int HCOD::sizeA() const
  {
    int s=0;
    for (size_t i=0; i<stages.size(); ++i)
      s+= stages[i]->sizeA();
    return s;
  }

  int HCOD::rank() const
  {
    int r=0;
    for (size_t i=0; i<stages.size(); ++i)
      r+= stages[i]->rank();
    return r;
  }

  void HCOD::setNameByOrder( const std::string root )
  {
    for (size_t i=0; i<stages.size(); ++i)
      {
	std::ostringstream os; os<<root<<i;
	stages[i]->name = os.str();
      }
  }

  /* --- DECOMPOSITION ------------------------------------------------------- */
  /* --- DECOMPOSITION ------------------------------------------------------- */
  /* --- DECOMPOSITION ------------------------------------------------------- */

  void HCOD::
  reset( void )
  {
    isReset=true;
    isInit=false;
    isSolutionCpt=false;

    solution.setZero(); Ytu.setZero();
    for( stage_iter_t iter = stages.begin();iter!=stages.end();++iter )
      {   (*iter)->reset();   }
  }

  void HCOD::
  initialize( void )
  {
    if(! isReset) reset(); // TODO: should it be automatically reset?
    assert( isReset&&(!isInit) ); isInit=true;

    /* Compute the initial COD for each stage. */
    unsigned int previousRank = 0;
    for( unsigned int i=0;i<stages.size();++i )
      {
	assert( stages[i]!=0 );
	sotDEBUG(5) <<" --- STAGE " <<i
		    << " --------------------------------- " << std::endl;
	previousRank
	  = stages[i]->computeInitialCOD(previousRank,soth::Stage::allRows(),Y);
      }
    Y.computeExplicitly(); // TODO: this should be done automatically on Y size.
  }
  void HCOD::
  update( const unsigned int & stageUp,const Stage::ConstraintRef & cst )
  {
    assert(isInit);
    GivensSequence Yup;
    unsigned int rankDef = stages[stageUp]->update(cst,Yup);
    for( unsigned int i=stageUp+1;i<stages.size();++i )
      {
	stages[i]->propagateUpdate(Yup,rankDef);
      }
    updateY(Yup);
  }
  void HCOD::
  update( stage_iter_t stageIter,const Stage::ConstraintRef & cst )
  {
    assert(isInit);
    sotDEBUG(5) << "Update " << (*stageIter)->name << std::endl;
    GivensSequence Yup;
    unsigned int rankDef = (*stageIter)->update(cst,Yup);
    for( ;stageIter!=stages.end();++stageIter )
      {
	(*stageIter)->propagateUpdate(Yup,rankDef);
      }
    updateY(Yup);
  }
  void HCOD::
  downdate( const unsigned int & stageDown, const unsigned int & rowDown )
  {
    assert(isInit);
    GivensSequence Ydown;
    bool propag=stages[stageDown]->downdate(rowDown,Ydown);
    for( unsigned int i=stageDown+1;i<stages.size();++i )
      {
     	propag = stages[i]->propagateDowndate(Ydown,propag);
      }
    updateY(Ydown);
  }
  void HCOD::
  downdate( stage_iter_t stageIter,const unsigned int & rowDown )
  {
    assert(isInit);
    sotDEBUG(5) << "Downdate " << (*stageIter)->name << std::endl;
    GivensSequence Ydown;
    bool propag=(*stageIter)->downdate(rowDown,Ydown);
    for( ;stageIter!=stages.end();++stageIter )
      {
     	propag = (*stageIter)->propagateDowndate(Ydown,propag);
      }
    updateY(Ydown);
  }



  void HCOD::
  updateY( const GivensSequence& Yup )
  {
    Y *= Yup;
    Yup.applyTransposeOnTheRight(Ytu);
    Yup.applyTransposeOnTheRight(Ytdu); /* TODO: this second multiplication could
					 * be avoided. */
  }


  /* --- COMPUTE ------------------------------------------------------------ */
  /* --- COMPUTE ------------------------------------------------------------ */
  /* --- COMPUTE ------------------------------------------------------------ */

  /* Assume the variable 'Ytu' contains Y^T*u.
   * The result is stored in the stages (getLagrangeMultipliers()).
   */
  void HCOD::
  computeLagrangeMultipliers()
  {
    assert( isSolutionCpt );
    stages.back()->computeRho(Ytu,rho,true);
    sotDEBUG(5) << "rho = " << (MATLAB)rho << std::endl;
    for( stage_riter_t iter=stages.rbegin()+1;iter!=stages.rend();++iter )
      {
	(*iter)->computeLagrangeMultipliers(rho);
      }
  }

  void HCOD::
  computeSolution(  bool compute_u )
  {
    assert(isInit);

    /* Initial solve. */
    Ytdu.setZero(); /* Ytdu.head(nullspace) only could be set to 0.
		     *  Does it make any diff? */
    for( unsigned int i=0;i<stages.size();++i )
      {
	stages[i]->computeSolution(Ytu,Ytdu,!isSolutionCpt);
      }
    sotDEBUG(5) << "Ytu = " << (MATLAB)Ytdu << std::endl;
    if( compute_u )
      {
	Y.multiply(Ytdu,du);
	sotDEBUG(5) << "u = " << (MATLAB)du << std::endl;
      }

    isSolutionCpt=true;
  }

  double HCOD::
  computeStepAndUpdate( void )
  {
    assert(isSolutionCpt);

    double tau = 1.0; Stage::ConstraintRef cst;
    stage_iter_t stageUp;
    for( stage_iter_t iter = stages.begin(); iter!=stages.end(); ++iter )
      {
	sotDEBUG(5) << "Check stage " << (*iter)->name << "." << std::endl;
	if(! (*iter)->checkBound( solution,du,cst,tau ) )
	  stageUp=iter;
      }
    if( tau<1 )
      {
	update(stageUp,cst);
     }
    return tau;
  }

  void HCOD::
  makeStep( double tau, bool compute_u )
  {
    Ytu += tau*Ytdu;
    if( compute_u ) { solution += tau*du; }
  }
  void HCOD::
  makeStep( bool compute_u )
  {
    Ytu += Ytdu;
    if( compute_u ) { solution += du; }
  }

  bool HCOD::
  searchAndDowndate( void )
  {
    double lambdamax = 0; unsigned int row;
    stage_iter_t stageDown;
    for( stage_iter_t iter = stages.begin(); iter!=stages.end(); ++iter )
      {
	if( (*iter)->maxLambda( lambdamax,row ) )
	  stageDown=iter;
      }
    if( lambdamax!=0 )
      {
	downdate(stageDown,row);
      }

  }


  /* --- ACTIVE SEARCH ------------------------------------------------------ */
  /* --- ACTIVE SEARCH ------------------------------------------------------ */
  /* --- ACTIVE SEARCH ------------------------------------------------------ */


    /* TODO:
       - Finish updateY

       - Compact the final active set (init aset is suppose to **TODO**
       be row-compact). / Assert this hypothesis.
       - Make all the necessary functions public, and externalize the algo. **TODO**

       - computation of the violation per stages **DONE**
       - computation of tau **DONE**
       - translation of cst_ref in triple<stage_ref,cst_ref,bound_ref> **NOT NECESSARY**
       - compute min lambda,w **DONE**
       - compute du from Ytdu and Ytu from u.  **STUPID**

       - Test with empty stages, full stages, rank 0 stages. **DONE**
       - Build a test with fixed-values matrices of all ranks. **DONE**

    */

  void HCOD::activeSearch( VectorXd & u )
  {
    /*
     * foreach stage: stage.initCOD(Ir_init)
     * u = 0
     * u0 = solve
     * do
     *   tau,cst_ref = max( violation(stages) )
     *   u += du*tau
     *   if( tau<1 )
     *     update(cst_ref); break;
     *
     *   lambda,w = computeLambda
     *   cst_ref,lmin = min( lambda,w )
     *   if lmin<0
     *     downdate( cst_ref )
     *
     */

    assert(VectorXi::LinSpaced(0,2,3)[0] == 0
            && VectorXi::LinSpaced(0,2,3)[1] == 1
            && VectorXi::LinSpaced(0,2,3)[2] == 2
            && "new version of Eigen might have change the "
	   "order of arguments in LinSpaced, please correct");

    initialize();
    bool endCondition = true;
    do
      {
	computeSolution();
	double tau = computeStepAndUpdate();
	if( tau<1 )
	  {
	    sotDEBUG(5) << "Update done, make step <1." << std::endl;
	    makeStep(tau);
	  }
	else
	  {
	    makeStep();
	    sotDEBUG(5) << "No update, make step <1." << std::endl;
	    computeLagrangeMultipliers();
	    if( searchAndDowndate() )
	      {
		sotDEBUG(5) << "Lagrange<0, downdate done." << std::endl;
	      }
	    else
	      {
		sotDEBUG(5) << "Lagrange>=0, no downdate." << std::endl;
		endCondition = false;
	      }
	  }
      } while(endCondition);

    u=solution;
  }



  /* --- TESTS -------------------------------------------------------------- */
  /* --- TESTS -------------------------------------------------------------- */
  /* --- TESTS -------------------------------------------------------------- */
  bool HCOD::
  testRecomposition( std::ostream* os )
  {
    bool res = true;
    for( unsigned int i=0;i<stages.size();++i )
      {
	bool sres=stages[i]->testRecomposition();
	if( os&&(!sres) ) *os << "Stage " <<i<<" is not properly recomposed."<<std::endl;
	res&=sres;
      }
    return res;
  }

  bool HCOD::
  testLagrangeMultipliers( std::ostream* os ) const
  {
    VectorXd verifL(sizeProblem); verifL.setZero();
    for( int i=0;i<stages.size();++i )
      {
	const Stage s = *stages[i];
	MatrixXd J_(s.nbConstraints(),sizeProblem);
	verifL += s.Jactive(J_).transpose()*s.getLagrangeMultipliers();
      }
    sotDEBUG(5) << "verif = " << (soth::MATLAB)verifL << std::endl;

    const double n = verifL.norm();
    const bool res = n<1e-6;
    if(os&&(!res)) (*os) << "TestLagrangian failed: norm is " << n << "."<<std::endl;
    return res;
  }

  void HCOD::
  show( std::ostream& os, bool check )
  {
    for( unsigned int i=0;i<stages.size();++i )
      {
	stages[i]->show(os,i+1,check);
      }

    MatrixXd Yex(sizeProblem,sizeProblem); Yex.setIdentity();
    Y.applyThisOnTheLeft(Yex);
    os<<"Y = " << (MATLAB)Yex << std::endl;
    if( isSolutionCpt )
      {
	os << "u = " << (MATLAB)solution << std::endl;
	os << "du = " << (MATLAB)du << std::endl;
      }
  }






}; // namespace soth

