#define SOTH_DEBUG
#define SOTH_DEBUG_MODE 15
#include "soth/debug.hpp"
#include "soth/COD.hpp"  // DEBUG
#include "soth/HCOD.hpp"
#include <boost/foreach.hpp>

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
    ,freezedStages(0)
    ,isReset(false),isInit(false),isSolutionCpt(false)
  {
    stages.reserve(nbStage);
    initialActiveSets.reserve(nbStage);
  }

  /* --- SETTERS/GETTERS ---------------------------------------------------- */
  /* --- SETTERS/GETTERS ---------------------------------------------------- */
  /* --- SETTERS/GETTERS ---------------------------------------------------- */

  void HCOD::
  pushBackStage( const MatrixXd & J, const VectorBound & bounds )
  {
    unsigned int s = stages.size();
    stages.resize( s+1 );
    stages[s] = stage_ptr_t(new soth::Stage( J,bounds,Y ));

    initialActiveSets.resize( s+1 );
    initialActiveSets[s].resize(0);
    isInit=false;
  }
  void HCOD::
  pushBackStage( const MatrixXd & J, const VectorBound & bounds,const VectorXi& Ir0 )
  {
    pushBackStage(J,bounds);
    setInitialActiveSet(Ir0,stages.size()-1);
  }
  void HCOD::
  pushBackStages( const std::vector<MatrixXd> & J,
		  const std::vector<VectorBound> & bounds )
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


  void HCOD::setDamping( const double & d )
  {
    for( stage_iter_t iter = stages.begin();iter!=stages.end();++iter )
      {	(*iter)->damping(d);      }
  }
  double HCOD::getMaxDamping()
  {
    double maxD=-1;
    for( stage_iter_t iter = stages.begin();iter!=stages.end();++iter )
      {
	double d = (*iter)->damping();
	assert( d>=0 );
	if( d>maxD ) maxD=d;
      }
    return maxD;
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

    Y.reset();
    solution.setZero(); Ytu.setZero();
    for( stage_iter_t iter = stages.begin();iter!=stages.end();++iter )
      {   (*iter)->reset();   }
  }

  void HCOD::
  initialize( void )
  {
    if(! isReset) reset(); // TODO: should it be automatically reset?
    assert( isReset&&(!isInit) );

    /* Compute the initial COD for each stage. */
    for( unsigned int i=0;i<stages.size();++i )
      {
	assert( stages[i]!=0 );
	sotDEBUG(5) <<" --- STAGE " <<i
		    << " ---------------------------------" << std::endl;
	stages[i]->computeInitialCOD(soth::Stage::allRows(),Y);
      }
    isReset=false; isInit=true;
  }
  void HCOD::
  update( const unsigned int & stageUp,const ConstraintRef & cst )
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
  update( stage_iter_t stageIter,const ConstraintRef & cst )
  {
    assert(isInit);
    sotDEBUG(5) << "Update " << (*stageIter)->name <<", "
		<< cst << std::endl;
    GivensSequence Yup;
    unsigned int rankDef = (*stageIter)->update(cst,Yup);
    for( ++stageIter;stageIter!=stages.end();++stageIter )
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
    sotDEBUG(5) << "Downdate " << (*stageIter)->name<<", row "<<rowDown
		<<" (cst="<<(*stageIter)->which(rowDown)<<")."<< std::endl;
    GivensSequence Ydown;
    bool propag=(*stageIter)->downdate(rowDown,Ydown);
    for( ++stageIter;stageIter!=stages.end();++stageIter )
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
  computeLagrangeMultipliers( const unsigned int & stageRef )
  {
    assert( isSolutionCpt );
    assert( stageRef<=stages.size() );

    stage_iter_t iter=stages.begin()+stageRef;

    if( iter==stages.end() )
      {   rho = - Ytu;     }
    else
      {
	(*iter)->computeRho(Ytu,rho,true); // This is Ytrho, not rho.
      }
    sotDEBUG(5) << "rho = " << (MATLAB)rho << std::endl;

    while ( iter!=stages.begin() )
      {
	iter--;
	(*iter)->computeLagrangeMultipliers(rho);
      }
  }

  void HCOD::
  computeSolution(  bool compute_u )
  {
    assert(isInit);

    /* Initial solve. */
    /* TODO: Ytdu.head(nullspace) only could be set to 0/Ytu.
     *  Does it make any diff? */
    if( isSolutionCpt ) Ytdu = -Ytu; else Ytdu.setZero();
    for( unsigned int i=0;i<stages.size();++i )
      {
	stages[i]->computeSolution(Ytu,Ytdu,!isSolutionCpt);
      }

    if( compute_u )
      {
	Y.multiply(Ytdu,du);
	sotDEBUG(5) << "u = " << (MATLAB)solution << std::endl;
	sotDEBUG(5) << "du = " << (MATLAB)du << std::endl;
      }

    sotDEBUG(5) << "Ytu = " << (MATLAB)Ytu << std::endl;
    sotDEBUG(5) << "Ytdu = " << (MATLAB)Ytdu << std::endl;
    isSolutionCpt=true;
  }
  void HCOD::
  damp( void )
  {
    assert(isInit);
    BOOST_FOREACH( stage_ptr_t sptr,stages )
      {	sptr->damp();      }
  }

  double HCOD::
  computeStepAndUpdate( void )
  {
    assert(isSolutionCpt);

    double tau = 1.0; ConstraintRef cst;
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
  double HCOD::
  computeStep( void )
  {
    assert(isSolutionCpt);

    double tau = 1.0; ConstraintRef cst;
    stage_iter_t stageUp;
    for( stage_iter_t iter = stages.begin(); iter!=stages.end(); ++iter )
      {
	sotDEBUG(5) << "Check stage " << (*iter)->name << "." << std::endl;
	if(! (*iter)->checkBound( solution,du,cst,tau ) )
	  stageUp=iter;
      }
    return tau;
  }

  void HCOD::
  makeStep( double tau, bool compute_u )
  {
    Ytu += tau*Ytdu;
    if( compute_u ) { solution += tau*du; }
    //makeStep(compute_u);
  }
  void HCOD::
  makeStep( bool compute_u )
  {
    Ytu += Ytdu;
    if( compute_u ) { solution += du; }
  }

  /* Return true iff the search is positive, ie if downdate was
   * needed and performed. */
  bool HCOD::
  searchAndDowndate( const unsigned int & stageRef )
  {
    assert( stageRef<=stages.size() );
    double lambdamax = Stage::EPSILON; unsigned int row;

    const stage_iter_t refend = stages.begin()+std::min(stages.size(),stageRef+1);
    stage_iter_t stageDown =  refend;
    for( stage_iter_t iter = stages.begin(); iter!=refend; ++iter )
      {
	if( (*iter)->maxLambda( solution,lambdamax,row ) )
	  stageDown=iter;
      }
    if( refend!=stageDown )
      {
	downdate(stageDown,row);
	return true;
      }
    return false;
  }

  bool HCOD::
  search( const unsigned int & stageRef )
  {
    assert( stageRef<=stages.size() );
    double lambdamax = Stage::EPSILON; unsigned int row;
    const stage_iter_t refend = stages.begin()+std::min(stages.size(),stageRef+1);
    stage_iter_t stageDown = refend;
    for( stage_iter_t iter = stages.begin(); iter!=refend; ++iter )
      {
	if( (*iter)->maxLambda( solution,lambdamax,row ) )
	  stageDown=iter;
      }
    return ( refend!=stageDown );
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
    sotDEBUGIN(15);
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

    assert(VectorXi::LinSpaced(3,0,2)[0] == 0
	   && VectorXi::LinSpaced(3,0,2)[1] == 1
	   && VectorXi::LinSpaced(3,0,2)[2] == 2
	   && "new version of Eigen might have change the "
	   "order of arguments in LinSpaced, please correct");

    initialize();
    Y.computeExplicitly(); // TODO: this should be done automatically on Y size.

    int iter = 0;
    unsigned int stageMinimal = 0;
    do
      {
	iter ++; sotDEBUG(5) << " --- *** \t" << iter << "\t***.---" << std::endl;

	if( sotDEBUG_ENABLE(15) )  show( sotDEBUGFLOW );
	assert( testRecomposition(&std::cerr) );
	//damp();
	computeSolution();
	assert( testSolution(&std::cerr) );

	double tau = computeStepAndUpdate();
	if( tau<1 )
	  {
	    sotDEBUG(5) << "Update done, make step <1." << std::endl;
	    makeStep(tau);
	  }
	else
	  {
	    sotDEBUG(5) << "No update, make step ==1." << std::endl;
	    makeStep();

	    for( ;stageMinimal<=stages.size();++stageMinimal )
	      {
		sotDEBUG(5) << "--- Started to examinate stage " << stageMinimal << std::endl;
		computeLagrangeMultipliers(stageMinimal);
		if( sotDEBUG_ENABLE(15) )  show( sotDEBUGFLOW );
		assert( testLagrangeMultipliers(stageMinimal,std::cerr) );

		if( searchAndDowndate(stageMinimal) )
		  {
		    sotDEBUG(5) << "Lagrange<0, downdate done." << std::endl;
		    break;
		  }

		for( unsigned int i=0;i<stageMinimal;++i )
		  stages[i]->freezeSlacks(false);
		if( stageMinimal<nbStages() )
		  stages[stageMinimal]->freezeSlacks(true);
	      }
	  }
    } while(stageMinimal<=nbStages());
    sotDEBUG(5) << "Lagrange>=0, no downdate, active search completed." << std::endl;

    u=solution;
    sotDEBUGOUT(15);
  }



  /* --- TESTS -------------------------------------------------------------- */
  /* --- TESTS -------------------------------------------------------------- */
  /* --- TESTS -------------------------------------------------------------- */
  bool HCOD::
  testRecomposition( std::ostream* os )
  {
    sotDEBUGPRIOR(+20);
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
  testSolution( std::ostream* os )
  {
    sotDEBUGPRIOR(+20);
    bool res = true;
    for( unsigned int i=0;i<stages.size();++i )
      {
	bool sres=stages[i]->testSolution( solution+du );
	if( os&&(!sres) ) *os << "Stage " <<i<<" has not been properly inverted."<<std::endl;
	res&=sres;
      }
    return res;



  }

  /* Compute sum(i=1:sr) Ji' li, with Jsr = I and lsr = u, and check
   * that the result is null. */
  bool HCOD::
  testLagrangeMultipliers( int unsigned stageRef,std::ostream* os ) const
  {
    assert( stageRef<=stages.size() );
    VectorXd verifL(sizeProblem);

    /* verifL = Jsr' lsr, with Jsr = I and lsr = u. */
    if( stageRef==stages.size() )
      { verifL = solution; stageRef -- ; }
    else
      verifL.setZero();

    /* verif += sum Ji' li. */
    for( unsigned int i=0;i<=stageRef;++i )
      {
	const Stage s = *stages[i];
	MatrixXd J_(s.nbConstraints(),sizeProblem);
	verifL += s.Jactive(J_).transpose()*s.getLagrangeMultipliers();
      }
    sotDEBUG(5) << "verif = " << (soth::MATLAB)verifL << std::endl;

    const double sumNorm = verifL.norm();
    const bool res = sumNorm<1e-6;
    if(os&&(!res))
      (*os) << "TestLagrangian failed: norm is " << sumNorm << "."<<std::endl;

    return res;
  }

  void HCOD::
  show( std::ostream& os, bool check )
  {
    sotDEBUGIN(15);
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
	os << "Ytu = " << (MATLAB)Ytu << std::endl;
	os << "Ytdu = " << (MATLAB)Ytdu << std::endl;
	assert( (solution-Y.matrixExplicit*Ytu).norm() < Stage::EPSILON );
	assert( (du-Y.matrixExplicit*Ytdu).norm() < Stage::EPSILON );
      }
    sotDEBUGOUT(15);
  }

  void HCOD::
  showActiveSet( std::ostream& os ) const
  {
    sotDEBUGIN(15);
    os << "{" << std::endl;
    for( unsigned int i=0;i<stages.size();++i )
      {
	os<< "    "; stages[i]->showActiveSet(os); os << std::endl;
      }
    os << "}" << std::endl;
    sotDEBUGOUT(15);
  }






}; // namespace soth

