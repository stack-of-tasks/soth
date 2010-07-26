/*
 *  Copyright
 *
 * Exhaustive tests of the HCOD, initial decompo, update and downdate.
 *
 */

#include "soth/HCOD.hpp"
#include "soth/debug.h"
#include "MatrixRnd.h"
#include <iomanip>


void generateDataSet( std::vector<Eigen::MatrixXd> &J,
		      std::vector<soth::bound_vector_t> &b,
		      const int NB_STAGE,
		      const int RANK[],
		      const int NR[],
		      const int NC )
{
  /* Initialize J and b. */
  J.resize(NB_STAGE);
  b.resize(NB_STAGE);

  unsigned int s = 0;
  for( int s=0;s<NB_STAGE;++s )
    {
      Eigen::MatrixXd Xhi( NR[ s],RANK[ s] );
      Eigen::MatrixXd Jfr( RANK[ s],NC );
      b[ s].resize(NR[ s]);

      soth::MatrixRnd::randomize( Xhi );
      soth::MatrixRnd::randomize( Jfr );
      J[s] = Xhi*Jfr;

      for( unsigned int i=0;i<NR[s];++i ) b[s][i] = (double)(i+1);
    }
}

void generateDeficientDataSet( std::vector<Eigen::MatrixXd> &J,
			       std::vector<soth::bound_vector_t> &b,
			       const int NB_STAGE,
			       const int RANKFREE[],
			       const int RANKLINKED[],
			       const int NR[],
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
      {
	Eigen::MatrixXd Xhifree( NR[s],RANKFREE[s] );
	Eigen::MatrixXd Jfr( RANKFREE[s],NC );
	soth::MatrixRnd::randomize( Xhifree );
	soth::MatrixRnd::randomize( Jfr );
	if( Xhifree.cols()>0 ) J[s] += Xhifree*Jfr;
      }
      if( (s>0)&&(RANKLINKED[s]>0) )
	{
	  Eigen::MatrixXd Xhilinked( NR[s],RANKLINKED[s] );
	  soth::MatrixRnd::randomize( Xhilinked );
	  Eigen::MatrixXd Alinked( RANKLINKED[s],NR[s-1] );
	  soth::MatrixRnd::randomize( Xhilinked );
	  J[s] += Xhilinked*Alinked*J[s-1];
	}

      for( unsigned int i=0;i<NR[s];++i ) b[s][i] = (double)(i+1);
    }
}

bool
clearIteralively( soth::HCOD& hcod )
{
  bool exitOk = true;
  for( int s=0;s<hcod.nbStages();++s )
    {
      while( hcod[s].sizeA()>0 )
	{
	  hcod.downdate(s,0);
	  exitOk&=hcod.testRecomposition(&std::cout);
	}
    }
  return exitOk;
}

int main (int argc, char** argv)
{
  sotDebugTrace::openFile();
  bool exitOk=true;

  Eigen::MatrixXd Massert(5,9);
  soth::MatrixRnd::randomize( Massert );
  assert( std::abs(Massert(1,2)-0.985007)<1e-5 );

  {
    /* All matrices full rank and independant, but the last one due to
     * lack of DOF.
     */
    const int NB_STAGE = 3;
    const int RANK[] = { 3,3,6 };
    const int NR[] = { 3,3,6 };
    const int NC = 10;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDataSet( J,b, NB_STAGE,RANK,NR,NC );

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();
    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );
  }

  {
    /* All matrices full rank, updating rows until the last rank saturates
     * du to the matrix size. Then remove all the lines from the first.
     */
    const int NB_STAGE = 3;
    const int RANK[] = { 3,3,6 };
    const int NR[] = { 3,3,6 };
    const int NC = 10;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDataSet( J,b,NB_STAGE,RANK,NR,NC );
    for( int i=2;i<NR[2];++i ) b[2][i] = std::make_pair(-i,i);

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    for( int i=2;i<NR[2];++i )
      { hcod.update( 2,std::make_pair(i,soth::Bound::BOUND_INF) ); }

    exitOk&=hcod.testRecomposition(&std::cout);
    //    if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    exitOk &= clearIteralively(hcod);
  }

  {
    /* All matrices rank def due to previous stages and by themselves.
     */
    const int NB_STAGE = 3;
    const int RANKFREE[] = { 2,2,3 };
    const int RANKLINKED[] = { 0,2,2 };
    const int NR[] = { 4,5,6 };
    const int NC = 12;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDeficientDataSet( J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC );
    //    for( unsigned int i=2;i<NR[2];++i ) b[2][i] = std::make_pair(-i,i);

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    exitOk &= clearIteralively(hcod);
  }


  exit( (exitOk?0:1) );

}



/* BUG LIST
 *
 *   - Insertion of a one-line matrix in the init: crash in the decompo QR.
 */
