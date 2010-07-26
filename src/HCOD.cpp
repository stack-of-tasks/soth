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
  {
    stages.reserve(nbStage);
    initialActiveSets.reserve(nbStage);
  }

  void HCOD::
  pushBackStage( const MatrixXd & J, const bound_vector_t & bounds )
  {
    unsigned int s = stages.size();
    stages.resize( s+1 );
    stages[s] = stage_ptr_t(new soth::Stage( J,bounds,Y ));

    initialActiveSets.resize( s+1 );
    initialActiveSets[s].resize(0);
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

  void HCOD::
  reset( void )
  {
    solution.setZero();
    throw "TODO";
  }

  void HCOD::
  initialize( void )
  {
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
    GivensSequence Yup;
    unsigned int rankDef = stages[stageUp]->update(cst,Yup);
    for( unsigned int i=stageUp+1;i<stages.size();++i )
      {
	stages[i]->propagateUpdate(Yup,rankDef);
      }
    updateY(Yup);
  }
  void HCOD::
  downdate( const unsigned int & stageDown, const unsigned int & rowDown )
  {
    GivensSequence Ydown;
    bool propag=stages[stageDown]->downdate(rowDown,Ydown);
    for( unsigned int i=stageDown+1;i<stages.size();++i )
      {
     	propag = stages[i]->propagateDowndate(Ydown,propag);
      }
    updateY(Ydown);
  }

  void HCOD::
  solve( void )
  {
    /* Compute the initial COD for each stage. */
    unsigned int previousRank = 0;
    for( unsigned int i=0;i<stages.size();++i )
      {
	assert( stages[i]!=0 );

	sotDEBUG(5) <<" --- STAGE " <<i<< " --------------------------------- " << std::endl;
	previousRank = stages[i]->computeInitialCOD(previousRank,soth::Stage::allRows(),Y);
      }

    /* Initial solve. */
    solution.setZero();
    for( unsigned int i=0;i<stages.size();++i )
      {
	stages[i]->solve(solution);
      }
    sotDEBUG(5) << "Ytu = " << (MATLAB)solution << std::endl;
    Y.applyThisOnVector( solution );
    sotDEBUG(5) << "u = " << (MATLAB)solution << std::endl;

    show(std::cout,true);
    Y.computeExplicitly();

    computeLambda();
    std::cout << "lambda = " << (MATLAB)lambda << std::endl;
    std::cout << "err = " << std::endl;
    for (int i=0; i<stages.size(); ++i)
      std::cout << "   " << (MATLAB)stages[i]->computeErr(solution) << std::endl;

    std::cout << " === DOWN ================================ " << std::endl;
    const unsigned int TO_DOWN = 0;
    const unsigned int ROW_DOWN = 0;
    downdate(TO_DOWN,ROW_DOWN);
    show(std::cout,true);

    std::cout << " === UP ================================ " << std::endl;
    const unsigned int TO_UP = 0;
    const unsigned int CSTR_UP = 5;
    update( TO_UP,std::make_pair(CSTR_UP,Bound::BOUND_INF) );
    show(std::cout,true);

  }

  void HCOD::
  updateY( const GivensSequence& Yup )
  {
    Y *= Yup;
    // TODO: update Ytu.
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

  //assume the variable 'solution' contains Y^T*u
  void HCOD::computeLambda()
  {
    int s = sizeA();
    lambda.resize(s);
    //last lambda
    int sn = stages.back()->sizeA();
    lambda.tail(sn) = stages.back()->computeErr(solution);
    s -= sn;
    VectorXd ro = stages.back()->computeRo(solution);
    assert(ro.size() == rank() - stages.back()->rank());
    int r = ro.size();

    for (int i=stages.size()-2; i>=0; --i)
    {
      int si = stages[i]->sizeA();
      s -= si;

      VectorBlock<VectorXd> l=lambda.segment(s, si);
      VectorBlock<VectorXd> roh = ro.head(r);
      stages[i]->computeLagrangeMultipliers(l,roh);
      r -= stages[i]->rank();
    }
  }


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

  void HCOD::
  show( std::ostream& os, bool check )
  {
    for( unsigned int i=0;i<stages.size();++i )
      {
	stages[i]->show(os,i,check);
      }

    MatrixXd Yex(sizeProblem,sizeProblem); Yex.setIdentity();
    Y.applyThisOnTheLeft(Yex);
    os<<"Y = " << (MATLAB)Yex << std::endl;
  }


  template< typename VectorGen >
  void HCOD::activeSearch( VectorGen & u )
  {
    assert(VectorXi::LinSpaced(0,2,3)[0] == 0
            && VectorXi::LinSpaced(0,2,3)[1] == 1
            && VectorXi::LinSpaced(0,2,3)[2] == 2
            && "new version of Eigen might have change the order of arguments in LinSpaced, please correct");
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

    /* TODO:
       - computation of the violation per stages
       - computation of tau
       - translation of cst_ref in triple<stage_ref,cst_ref,bound_ref>
       - compute du from Ytdu and Ytu from u.
       - compute min lambda,w
       - Make all the necessary functions public, and externalize the algo.

       - Test with empty stages, full stages, rank 0 stages.
       - Build a test with fixed-values matrices of all ranks.

       - Compact the final active set (init aset is suppose to
       be row-compact). / Assert this hypothesis.
    */

  }





}; // namespace soth

