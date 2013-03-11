/* -------------------------------------------------------------------------- *
 * 
 * Simple test of the HCOD class with assertion.
 * 
 * -------------------------------------------------------------------------- */

#include "soth/debug.hpp"
#include "soth/HCOD.hpp"
#include "MatrixRnd.hpp"
#include "RandomGenerator.hpp"

int main (int argc, char** argv)
{
  soth::sotDebugTrace::openFile();
  const int NB_STAGE = 3;
  const int RANKFREE[]   = { 3, 4, 3,     5, 3 };
  const int RANKLINKED[] = { 2, 2, 1,     5, 3 };
  const int NR[]         = { 5, 4, 5,     5, 8 };
  const int NC = 12;

  /* Initialize J and b. */
  std::vector<Eigen::MatrixXd> J(NB_STAGE);
  std::vector<soth::VectorBound> b(NB_STAGE);
  soth::generateDeficientDataSet(J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC);
  b[0][1] = std::make_pair(-0.1,2.37);
  for( int i=0;i<NB_STAGE;++i )
    {
      sotDEBUG(0) << "J"<<i+1<<" = " << (soth::MATLAB)J[i] << std::endl;
      sotDEBUG(0) << "e"<<i+1<< " = " << b[i] << ";"<<std::endl;
    }
  assert( std::abs(J[0](0,0)-(-1.1149))<1e-5 );

  /* SOTH structure construction. */
  soth::HCOD hcod(NC,NB_STAGE);
  for( int i=0;i<NB_STAGE;++i )
    {
      hcod.pushBackStage( J[i],b[i] );
      assert(NR[i]>0);
    }
  hcod.setInitialActiveSet();
  hcod.setNameByOrder("stage_");

  hcod.initialize();
  hcod.computeSolution(true);
  assert(  (hcod.rank()==10)   && (hcod[0].rank()==3)
	 &&(hcod[1].rank()==4) && (hcod[2].rank()==3)  );

  if( soth::sotDEBUGFLOW.outputbuffer.good() )
    hcod.show( soth::sotDEBUGFLOW.outputbuffer );
  double tau = hcod.computeStepAndUpdate();
  hcod.makeStep(tau);
  assert((std::abs(tau-1.)<=10*soth::Stage::EPSILON)&&"Check bound test failed.");

  hcod.computeLagrangeMultipliers(hcod.nbStages());
  bool testL = hcod.testLagrangeMultipliers(hcod.nbStages(),std::cout);
  sotDEBUG(5) << "Test multipliers: " << ((testL)?"Passed!":"Failed...") << std::endl;
  assert(testL&&"Lagrange Multipliers test failed.");
}
