/*
 *  Copyright
 */
//#define SOTH_DEBUG
//#define SOTH_DEBUG_MODE 45
#include "soth/debug.hpp"
#include "soth/HCOD.hpp"
#include "soth/BasicStage.hpp"
#include "MatrixRnd.hpp"
#include <sys/time.h>
#include <Eigen/SVD>
#include "RandomGenerator.hpp"

using namespace soth;
using std::endl;

/* -------------------------------------------------------------------------- */
namespace soth
{
  struct Now
  {
    long int sec,usec;
    Now(void)
    {
      struct timeval tref;
      gettimeofday(&tref,NULL);
      sec=tref.tv_sec; usec=tref.tv_usec;
    }
  };
  double operator-( const Now& t1, const Now& t0 )
  {
    return (t1.sec-t0.sec)*1000.0 + (t1.usec-t0.usec+0.0)/1000.0;
  }
  std::ostream& operator<< (std::ostream& os,const Now& now )
  { return os << now.sec <<"' " << now.usec;  }

};



/* -------------------------------------------------------------------------- */

void testBasicStage()
{
  MatrixXd m1(5,4);
  Map<MatrixXd> map1(m1.data(), m1.size(), 1);
  map1 = VectorXd::LinSpaced(m1.size(), 0, m1.size()-1);
  std::cout << "m1 = " << m1 << endl;

  VectorBound b1(5);
  b1[0] = 1.6;
  b1[1] = std::make_pair(-0.1,0.2);
  std::cout << "b1 = " <<b1 << endl;

  soth::BaseY Y(5);
  soth::BasicStage st( 5,4,m1.data(),b1.data(),Y );
  st.set( m1.data(),b1.data() );

  std::cout << "J=" << st.getJ() << endl;
  std::cout << "b=" << st.getBound() << endl;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
int main (int argc, char** argv)
{
# ifndef NDEBUG
  sotDebugTrace::openFile();
#endif

  unsigned int NB_STAGE,NC;
  std::vector<unsigned int> NR,RANKLINKED,RANKFREE;
  std::vector<Eigen::MatrixXd> J;
  std::vector<soth::bound_vector_t> b;

  if( (argc==3)&& std::string(argv[1])=="-file")
    {
      readProblemFromFile( argv[2],J,b,NB_STAGE,NR,NC);
    }
  else
    {
      /* Initialize the seed. */
      struct timeval tv;
      gettimeofday(&tv,NULL);
      int seed = tv.tv_usec % 7919; //= 7594;
      if( argc == 2 )
	{  seed = atoi(argv[1]);  }
      std::cout << "seed = " << seed << std::endl;
      soth::Random::setSeed(seed);

      /* Decide the size of the problem. */
      generateRandomProfile(NB_STAGE,RANKFREE,RANKLINKED,NR,NC);
      /* Initialize J and b. */
      generateDeficientDataSet(J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC);
    }

  std::cout << "NB_STAGE=" << NB_STAGE <<",  NC=" << NC << endl;
  std::cout << endl;
  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      std::cout << "J"<<i+1<<" = " << (soth::MATLAB)J[i] << std::endl;
      std::cout << "e"<<i+1<< " = " << b[i] << ";"<<std::endl;
    }

  //if( seed==317)  assert( std::abs(J[0](0,0)-(-1.1149))<1e-5 );

  /* SOTH structure construction. */
  soth::HCOD hcod(NC,NB_STAGE);
  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      hcod.pushBackStage( J[i],b[i] );
    }
  hcod.setNameByOrder("stage_");



  testBasicStage();
}
