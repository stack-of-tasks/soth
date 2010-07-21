#include "soth/HCOD.hpp"



namespace soth
{


  HCOD::
  HCOD( unsigned int sizeProblem, unsigned int nbStage )
  :
    Y(sizeProblem)
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

	std::cout <<" --- STAGE " <<i<< " --------------------------------- " << std::endl;
	previousRank = stages[i]->computeInitialCOD(previousRank,soth::Stage::allRows());
      }

    for( unsigned int i=0;i<stages.size();++i )
      {
	Eigen::MatrixXd Jrec; stages[i]->recompose(Jrec);
	std::cout << "Jrec" <<i<<" = " << (soth::MATLAB)Jrec << std::endl;
      }


  }





}; // namespace soth

