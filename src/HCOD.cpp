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
    std::cout << "lambda = " << std::endl << lambda << std::endl;
    std::cout << "err = " << std::endl;
    for (int i=0; i<stages.size(); ++i)
      std::cout << stages[i]->computeErr(solution) << std::endl;

    std::cout << " === DOWN ================================ " << std::endl;
    const unsigned int TO_DOWN = 0;
    const unsigned int ROW_DOWN = 1;
    GivensSequence Ydown;

    bool propag=stages[TO_DOWN]->downdate(ROW_DOWN,Ydown);

    for( unsigned int i=TO_DOWN+1;i<stages.size();++i )
      {
     	propag = stages[i]->propagateDowndate(Ydown,propag);
      }
    updateY(Ydown);

    show(std::cout,true);

    std::cout << " === UP ================================ " << std::endl;
    const unsigned int TO_UP = 0;
    const unsigned int CSTR_UP = 5;
    GivensSequence Yup;
    unsigned int rankDef
      = stages[TO_UP]->update( std::make_pair(CSTR_UP,Bound::BOUND_TWIN),Yup);
    // TODO: propagate.
    for( unsigned int i=TO_UP+1;i<stages.size();++i )
      {
     	//TODO
	stages[i]->propagateUpdate(Ydown,rankDef);
      }
    updateY(Yup);

    show(std::cout,true);

    /* --- RECOMPOSE FOR TEST --- */
    // for( unsigned int i=0;i<stages.size();++i )
    //   {
    // 	Eigen::MatrixXd Jrec; stages[i]->recompose(Jrec);
    // 	sotDEBUG(5) << "Jrec" <<i<<" = " << (soth::MATLAB)Jrec << std::endl;
    //   }


  }

  void HCOD::
  updateY( const GivensSequence& Yup )
  {
    Y *= Yup;

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



}; // namespace soth

