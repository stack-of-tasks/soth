/*
 *  Copyright
 */
//#define SOTH_DEBUG
//#define SOTH_DEBUG_MODE 45
#include "soth/debug.h"
#include "soth/HCOD.hpp"
#include "soth/debug.h"
#include "MatrixRnd.h"
#include <sys/time.h>
#include <Eigen/SVD>

namespace Eigen
{
  #include "soth/DestructiveColPivQR.h"
}

using namespace soth;
using std::endl;

void generateDeficientDataSet( std::vector<Eigen::MatrixXd> &J,
			       std::vector<soth::bound_vector_t> &b,
			       const int NB_STAGE,
			       const VectorXd & RANKFREE,
			       const VectorXd & RANKLINKED,
			       const VectorXd & NR,
			       const int NC )
{
  /* Initialize J and b. */
  J.resize(NB_STAGE);
  b.resize(NB_STAGE);

  unsigned int s = 0;
  for( int s=0;s<NB_STAGE;++s )
    {
      b[ s].resize(NR[ s]);

      assert( (RANKFREE[s]>0)||(RANKLINKED[s]>0) );

      J[s].resize( NR[s],NC ); J[s].setZero();
      if( RANKFREE[s]>0 )
	{
	  Eigen::MatrixXd Xhifree( NR[s],RANKFREE[s] );
	  Eigen::MatrixXd Jfr( RANKFREE[s],NC );
	  soth::MatrixRnd::randomize( Xhifree );
	  soth::MatrixRnd::randomize( Jfr );
	  if( Xhifree.cols()>0 ) J[s] += Xhifree*Jfr;
	}
      if( RANKLINKED[s]>0 )
	{
	  Eigen::MatrixXd Xhilinked( NR[s],RANKLINKED[s] );
	  soth::MatrixRnd::randomize( Xhilinked );
	  for( int sb=0;sb<s;++sb )
	  {
	    Eigen::MatrixXd Alinked( RANKLINKED[s],NR[sb] );
	    soth::MatrixRnd::randomize( Alinked );
	    J[s] += Xhilinked*Alinked*J[sb];
	  }
	}

      for( unsigned int i=0;i<NR[s];++i ) b[s][i] = (double)(i+1);
    }
}


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


  // std::ostream & chronout = std::cout;

  // class chrono_context_ptr_t;
  // struct Chrono
  // {
  //   chrono_context_ptr_t context;
  //   std::string definition;
  //   Chrono( chrono_context_ptr_t context,std::string definition )
  //     : context(context),definition(definition)
  //   {}
  //   ~Chrono()
  //   {
  //     chronout << definition << ": " << (Now()-context->t0) << std::endl;
  //   }
  // };
  // boost::smart_ptr<Chrono> chrono_ptr_t;


  // struct ChronoContext
  // {
  //   time_t t0;
  //   std::list<chrono_ptr_t> chronoList;
  // };
  // std::map< std::string,ChronoContext >

};



/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
int main (int argc, char** argv)
{
  {
    struct timeval tv;
    gettimeofday(&tv,NULL);

    int seed = tv.tv_usec % 7919;
    if( argc == 2 )
      {  seed = atoi(argv[1]);  }
    std::cout << "seed = " << seed << std::endl;
    soth::Random::setSeed(seed);
  }
#ifdef SOTH_DEBUG
  sotDebugTrace::openFile();
#endif

  /* Decide the size of the problem. */
  // int NB_STAGE = 6,NC = 46;
  // Eigen::VectorXd NR(NB_STAGE),RANKLINKED(NB_STAGE),RANKFREE(NB_STAGE);
  // NR         << 6, 3, 6, 6, 6, 6;
  // RANKFREE   << 6, 3, 4, 4, 5, 6;
  // RANKLINKED << 0, 0, 2, 2, 1, 0;
  int NB_STAGE = 4,NC = 30;
  Eigen::VectorXd NR(NB_STAGE),RANKLINKED(NB_STAGE),RANKFREE(NB_STAGE);
  NR         << 6, 3, 6, 6;
  RANKFREE   << 6, 3, 6, 6;
  //  RANKFREE   << 6, 3, 4, 4;
  RANKLINKED << 0, 0, 2, 2;

  /* Initialize J and b. */
  std::vector<Eigen::MatrixXd> J(NB_STAGE);
  std::vector<soth::bound_vector_t> b(NB_STAGE);
  generateDeficientDataSet(J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC);
  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      //std::cout << "J"<<i+1<<" = " << (soth::MATLAB)J[i] << std::endl;
      //std::cout << "e"<<i+1<< " = " << b[i] << ";"<<std::endl;
    }
  //if( seed==317)  assert( std::abs(J[0](0,0)-(-1.1149))<1e-5 );

  /* SOTH structure construction. */
  soth::HCOD hcod(NC,NB_STAGE);
  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      hcod.pushBackStage( J[i],b[i] );
    }
  hcod.setNameByOrder("stage_");

  VectorXd solution;

  Now t0;
  struct timeval tv0,tv1;

  for( unsigned int i=0;i<1000;++i )
    {   hcod.reset(); hcod.initialize(); }

  Now t1;
  for( unsigned int i=0;i<1000;++i )
    hcod.Y.computeExplicitly();

  Now t2;
  for( unsigned int i=0;i<1000;++i )
    hcod.computeSolution(true);

  Now t3;
  for( unsigned int i=0;i<1000;++i )
    hcod.reset();

  Now t4;
  for( unsigned int i=0;i<1000;++i )
    double tau = hcod.computeStep();

  Now t5;
  for( unsigned int i=0;i<1000;++i )
    hcod.computeLagrangeMultipliers();

  Now t6;
  for( unsigned int i=0;i<1000;++i )
    bool down = hcod.search();

  Now t7;
  for( unsigned int i=0;i<1000;++i )
    hcod.activeSearch( solution );

  Now tf;
  std::cout << "HCOD = " <<(t1-t0) <<"us"<< endl;
  std::cout << "Y = " << (t2-t1) <<"us"<< endl;
  std::cout << "Inv = " << (t3-t2) <<"us"<< endl;
  std::cout << "reset = " << (t4-t3) <<"us"<< endl;
  std::cout << "check = " << (t5-t4) <<"us"<< endl;
  std::cout << "lagrange = " << (t6-t5) <<"us"<< endl;
  std::cout << "l<0 = " << (t7-t6) <<"us"<< endl;
  std::cout << "activesearch = " << (tf-t7) <<"us"<< endl;

  if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

  std::cout << " --- UNITCHRONO --------------------------------- " << endl;

  { // simple QR test
    const int Nctest = NC;
    const int Nrtest = 21;
    MatrixXd A = MatrixXd::Random(Nrtest,Nctest);
    MatrixXd Y(Nctest,Nctest);

    Now t3;
    for( int i=0;i<1000;++i )
      Eigen::DestructiveColPivQR<MatrixXd,MatrixXd> mQR(A,Y,1e-6);
    Now tf;
    std::cout << "QR " << Nctest << "x" << Nrtest << " = " << tf-t3 <<"us"<< endl;
  }

  { // simple Linv test
    const int Nctest = NC;
    const int Nrtest = 21;
    MatrixXd A = MatrixXd::Random(Nrtest,Nctest);
    VectorXd b(Nrtest);

    Now t3;
    for( int i=0;i<1000;++i )
      solveInPlaceWithLowerTriangular( A.leftCols(Nrtest),b );
    Now tf;
    std::cout << "Linv " << Nrtest << "x" << Nrtest << " = " << tf-t3 <<"us"<< endl;
  }


  { // Double test with matrix multiplication, just to be sure.
    const int Nctest = NC;
    const int Nrtest = 21;
    MatrixXd A = MatrixXd::Random(Nrtest,Nctest);
    MatrixXd B = MatrixXd::Random(Nctest,Nctest);

    MatrixXd C(Nrtest,Nctest);
    Now t3;
    for( int i=0;i<1000;++i ) C=A*B;
    Now tf;
    std::cout << Nctest << "x" << Nrtest << " = " << tf-t3 <<"us"<< endl;
    std::cout << C(0,0);
  }

}
