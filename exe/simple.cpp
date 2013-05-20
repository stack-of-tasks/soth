/* -------------------------------------------------------------------------- *
 * 
 * IJRR 2013 tests: compare the cost of the HQP with the QP cascade of Kanoun
 * and DeLassa.
 * 
 * -------------------------------------------------------------------------- */
#define SOTH_DEBUG
#define SOTH_DEBUG_MODE 45
#include "soth/debug.hpp"
#include "soth/HCOD.hpp"
#include "soth/Random.hpp"
#include "COD.hpp"
#include "RandomGenerator.hpp"

#include <boost/assign/std/vector.hpp> // for 'operator+=()'
using namespace boost::assign; // bring 'operator+=()' into scope

#ifndef WIN32
#include <sys/time.h>
#endif // WIN32

#include <iostream>
#include <sstream>
#include "gettimeofday.hpp"

using namespace soth;
using std::endl;
using std::cout;
using std::cerr;

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

void ehqp( std::vector<Eigen::MatrixXd> J,
	  std::vector<soth::VectorBound> b,
	  unsigned int NB_STAGE,
	  unsigned int NC )
{
  HCOD hsolver(NC,NB_STAGE); VectorXd solution(NC);
  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      hsolver.pushBackStage(J[i],b[i]);
      hsolver.setNameByOrder("level");
    }
  hsolver.setDamping(0.0);
  hsolver.setInitialActiveSet();

  hsolver.initialize();
  hsolver.Y.computeExplicitly();
  hsolver.computeSolution();
}
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

void ihqp( std::vector<Eigen::MatrixXd> J,
	  std::vector<soth::VectorBound> b,
	  unsigned int NB_STAGE,
	  unsigned int NC )
{
  HCOD hsolver(NC,NB_STAGE); VectorXd solution(NC);
  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      hsolver.pushBackStage(J[i],b[i]);
      hsolver.setNameByOrder("level");
    }
  hsolver.setDamping(0.0);
  hsolver.setInitialActiveSet();
  hsolver.activeSearch(solution);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */


int main ()
{
  unsigned int NB_STAGE,NC;
  std::vector<unsigned int> NR,RANKLINKED,RANKFREE;
  std::vector<Eigen::MatrixXd> J;
  std::vector<soth::VectorBound> b;

  soth::Random::setSeed(704819);
  const int size = 40; 

  generateFixedSizeRandomProfile(size,
				 1,0.8,1,NB_STAGE,RANKFREE,RANKLINKED,NR,NC);
  generateDeficientDataSet(J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC);

  struct timeval t0,t1;
  double time;

  gettimeofday(&t0,NULL);
  ehqp(J,b,NB_STAGE,NC);
  gettimeofday(&t1,NULL);
  time = (t1.tv_sec-t0.tv_sec)+(t1.tv_usec-t0.tv_usec)/1.0e6;
  cout << "ehqp = " << time*1000 << " ms " << std::endl;

  gettimeofday(&t0,NULL);
  ihqp(J,b,NB_STAGE,NC);
  gettimeofday(&t1,NULL);
  time = (t1.tv_sec-t0.tv_sec)+(t1.tv_usec-t0.tv_usec)/1.0e6;
  cout << "ihqp = " << time*1000 << " ms " << std::endl;

}
