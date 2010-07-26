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
  const int executeAll = 1;

  Eigen::MatrixXd Massert(5,9);
  soth::MatrixRnd::randomize( Massert );
  assert( std::abs(Massert(1,2)-0.985007)<1e-5 );

  if(executeAll){
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
    assert(exitOk);
  }

  if(executeAll){
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
    for( int i=2;i<NR[2];++i ) b[2][i] = std::make_pair(-i-1,i+1);

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    for( int i=2;i<NR[2];++i )
      { hcod.update( 2,std::make_pair(i,soth::Bound::BOUND_INF) ); }

    exitOk&=hcod.testRecomposition(&std::cout);
    //    if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }

  /* --- UPDATE TESTS ------------------------------------------------------- */
  if(executeAll){
    /* Insertion of a full rank line, increasing the total rank, on stages first, middle, last.
     */
    const int NB_STAGE = 3;
    const int RANK[] = { 3,3,3 };
    const int NR[] = { 3,3,3 };
    const int NC = 10;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDataSet( J,b,NB_STAGE,RANK,NR,NC );
    for( int s=0;s<NB_STAGE;++s ) b[s][2] = std::make_pair(-3,3);

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    int rank=hcod.rank();
    for( int i=0;i<NB_STAGE;++i )
      {
	hcod.update( i,std::make_pair(2,soth::Bound::BOUND_INF) );
	exitOk&=hcod.testRecomposition(&std::cout);
	assert( hcod.rank()==++rank );
     }

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }

  if(executeAll){
    /* Insertion of a rank-def line, preserving the total rank, on stages first, middle, last.
     */
    const int NB_STAGE = 3;
    const int RANK[] = { 2,2,2 };
    const int NR[] = { 3,3,3 };
    const int NC = 10;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDataSet( J,b,NB_STAGE,RANK,NR,NC );
    for( int s=0;s<NB_STAGE;++s ) b[s][2] = std::make_pair(-3,3);

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    int rank=hcod.rank();
    for( int i=0;i<NB_STAGE;++i )
      {
	const int rankStage = hcod[i].rank();
	hcod.update( i,std::make_pair(2,soth::Bound::BOUND_INF) );
	exitOk&=hcod.testRecomposition(&std::cout);
	assert( hcod.rank()==rank );
  	assert( rankStage == hcod[i].rank() );
    }

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }

  if(executeAll){
    /* Insertion of a full rank line, increasing stage rank but not the total rank,
     * on stages first, and middle (not last, this case is not possible).
     */
    const int NB_STAGE = 3;
    const int RANK[] = { 3,3,3 };
    const int NR[] = { 3,3,3 };
    const int NC = 10;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDataSet( J,b,NB_STAGE,RANK,NR,NC );
    for( int s=0;s<NB_STAGE-1;++s )
      {
	const int LAST = NR[s]-1;
	b[s][LAST] = std::make_pair(-LAST-1,LAST-1);
      }
    { // Link last line of 0 to first 2 lines of 1.
      const int s=0,LAST = NR[s]-1;
      Eigen::MatrixXd Xhi(1,NR[s+1]-1); soth::MatrixRnd::randomize(Xhi);
      J[s].row(LAST) = Xhi*J[s+1].topRows(NR[s+1]-1);
    }
    { // Link last line of 1 to 2.
      const int s=1,LAST = NR[s]-1;
      Eigen::MatrixXd Xhi(1,NR[s+1]); soth::MatrixRnd::randomize(Xhi);
      J[s].row(LAST) = Xhi*J[s+1];
    }

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    int rank=hcod.rank();
    for( int i=0;i<NB_STAGE-1;++i )
      {
	const int rankStage = hcod[i].rank();
	hcod.update( i,std::make_pair(2,soth::Bound::BOUND_INF) );
	exitOk&=hcod.testRecomposition(&std::cout);
	assert( hcod.rank()==rank );
	assert( rankStage+1 == hcod[i].rank() );
     }

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }

  /* --- Overshoot --- */
  if(executeAll){
    /* Insertion at the last stage when rank==nc (direct overshoot possible).
     */
    const int NB_STAGE = 5;
    const int RANK[] = { 3,3,3,3,5 };
    const int NR[] = { 4,4,4,4,5 };
    const int NC = 15;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDataSet( J,b,NB_STAGE,RANK,NR,NC );

    const int LS = NB_STAGE-1, LR=NR[LS]-1;
    b[LS][LR] = std::make_pair(-LR-1,LR+1);

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    const int rank=hcod.rank();
    const int rankStage = hcod[LS].rank();
    assert(rank==NC);
    hcod.update( LS,std::make_pair(LR,soth::Bound::BOUND_INF) );
    exitOk&=hcod.testRecomposition(&std::cout);
    assert( hcod.rank()==rank );
    assert( rankStage == hcod[LS].rank() );

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }

  if(executeAll){
    /* Insertion at the last stage when rank==nc (indirect overshoot possible).
     */
    const int NB_STAGE = 5;
    const int RANK[] = { 3,3,3,3,5 };
    const int NR[] = { 4,4,4,4,5 };
    const int NC = 15;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDataSet( J,b,NB_STAGE,RANK,NR,NC );

    const int LS = NB_STAGE-3, LR=NR[LS]-1;
    b[LS][LR] = std::make_pair(-LR-1,LR+1);

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    const int rank=hcod.rank();
    const int rankStage = hcod[LS].rank();
    assert(rank==NC);
    hcod.update( LS,std::make_pair(LR,soth::Bound::BOUND_INF) );
    exitOk&=hcod.testRecomposition(&std::cout);
    assert( hcod.rank()==rank );
    assert( rankStage == hcod[LS].rank() );

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }

  /* --- Closure --- */
  if(executeAll){
    /* Insertion of a full rank line that causes the closure of a stage above.
     */
    const int NB_STAGE = 3;
    const int RANKFREE[] = { 6,0,0 };
    const int RANKLINKED[] = { 0,2,3 };
    const int NR[] = { 6,2,3 };
    const int NC = 12;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDeficientDataSet( J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC );
    for( unsigned int i=2;i<NR[0];++i )
      b[0][i] = std::make_pair(-i-1,i+1);

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );
    assert( hcod.rank()==6 );       assert( hcod[0].rank()==2 );
    assert( hcod[1].rank()==2 );    assert( hcod[2].rank()==2 );

    hcod.update( 0,std::make_pair(2,soth::Bound::BOUND_INF) );
    exitOk&=hcod.testRecomposition(&std::cout);
    assert( hcod.rank()==6 );       assert( hcod[0].rank()==3 );
    assert( hcod[1].rank()==2 );    assert( hcod[2].rank()==1 );

    hcod.update( 0,std::make_pair(3,soth::Bound::BOUND_INF) );
    exitOk&=hcod.testRecomposition(&std::cout);
    assert( hcod.rank()==6 );       assert( hcod[0].rank()==4 );
    assert( hcod[1].rank()==2 );    assert( hcod[2].rank()==0 );

    // hcod.update( 0,std::make_pair(4,soth::Bound::BOUND_INF) );
    // exitOk&=hcod.testRecomposition(&std::cout);
    // assert( hcod.rank()==6 );       assert( hcod[0].rank()==5 );
    // assert( hcod[1].rank()==1 );    assert( hcod[2].rank()==0 );

    // TODO: Encore un bug a corriger ici
    // hcod.update( 0,std::make_pair(5,soth::Bound::BOUND_INF) );
    // exitOk&=hcod.testRecomposition(&std::cout);
    // assert( hcod.rank()==6 );       assert( hcod[0].rank()==6 );
    // assert( hcod[1].rank()==0 );    assert( hcod[2].rank()==0 );

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }


  // TODO: closing of a stage due to rank loss.

  /* --- DOWNDATE TESTS ----------------------------------------------------- */
  if(executeAll){
    /* Removal of a full rank line, with decrease of the total rank.
     */
    const int NB_STAGE = 3;
    const int RANK[] = { 3,3,3 };
    const int NR[] = { 3,3,3 };
    const int NC = 10;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDataSet( J,b,NB_STAGE,RANK,NR,NC );

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    int rank=hcod.rank();
    for( int i=0;i<NB_STAGE;++i )
      {
	const int rankStage = hcod[i].rank();
 	hcod.downdate( i,1 );
	exitOk&=hcod.testRecomposition(&std::cout);
	if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );
	assert( hcod.rank()==--rank );
	assert( rankStage-1 == hcod[i].rank() );
     }

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }

  if(executeAll){
    /* Removal of a rank def line, with no decrease of neither the stage nor the total rank.
     */
    const int NB_STAGE = 3;
    const int RANK[] = { 2,2,2 };
    const int NR[] = { 3,3,3 };
    const int NC = 10;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDataSet( J,b,NB_STAGE,RANK,NR,NC );

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    int rank=hcod.rank();
    for( int i=0;i<NB_STAGE;++i )
      {
	const int rankStage = hcod[i].rank();
 	hcod.downdate( i,1 );
	exitOk&=hcod.testRecomposition(&std::cout);
	if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );
	assert( hcod.rank()==rank );
	assert( rankStage == hcod[i].rank() );
     }

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }

  if(executeAll){
    /* Insertion of a full rank line, decreasing stage rank but not the total rank,
     * on stages first, and middle (not last, this case is not possible).
     */
    const int NB_STAGE = 3;
    const int RANK[] = { 3,3,3 };
    const int NR[] = { 3,3,3 };
    const int NC = 10;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDataSet( J,b,NB_STAGE,RANK,NR,NC );
    { // Link last line of 0 to first 2 lines of 1.
      const int s=0,LAST = NR[s]-1;
      Eigen::MatrixXd Xhi(1,NR[s+1]-1); soth::MatrixRnd::randomize(Xhi);
      J[s].row(LAST) = Xhi*J[s+1].topRows(NR[s+1]-1);
    }
    { // Link last line of 1 to 2.
      const int s=1,LAST = NR[s]-1;
      Eigen::MatrixXd Xhi(1,NR[s+1]); soth::MatrixRnd::randomize(Xhi);
      J[s].row(LAST) = Xhi*J[s+1];
    }

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );
    assert( hcod[0].rank()==3 );
    assert( hcod[1].rank()==2 );
    assert( hcod[2].rank()==2 );
    assert( hcod[0].gete()(hcod[0].where(2),0) == 3 );

    int rank=hcod.rank();
    for( int i=0;i<NB_STAGE-1;++i )
      {
	const int rankStage = hcod[i].rank();
	hcod.downdate( i,hcod[i].where(NR[i]-1) );
	exitOk&=hcod.testRecomposition(&std::cout);
	assert( hcod.rank()==rank );
	assert( rankStage-1 == hcod[i].rank() );
      }

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }

  /* --- RANK DEFICIENCY ---------------------------------------------------- */
  /* --- RANK DEFICIENCY ---------------------------------------------------- */
  /* --- RANK DEFICIENCY ---------------------------------------------------- */
  exit(0);

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
    //    for( unsigned int i=2;i<NR[2];++i ) b[2][i] = std::make_pair(-i-1,i+1);

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }

  {
    /* Same as before, but iterative construction ---> TODO
     */
    const int NB_STAGE = 3;
    const int RANKFREE[] = { 2,2,3 };
    const int RANKLINKED[] = { 0,2,2 };
    const int NR[] = { 4,5,6 };
    const int NC = 12;

    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDeficientDataSet( J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC );
    for( unsigned int s=0;s<NB_STAGE;++s )
      for( unsigned int i=2;i<NR[s];++i )
	b[s][i] = std::make_pair(-i-1,i+1);

    soth::HCOD hcod(NC,NB_STAGE);
    hcod.pushBackStages( J,b );

    hcod.initialize();

    exitOk&=hcod.testRecomposition(&std::cout);
    //if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );

    exitOk &= clearIteralively(hcod);
    assert(exitOk);
  }


  exit( (exitOk?0:1) );

}



/* BUG LIST
 *
 *   - Insertion of a one-line matrix in the init: crash in the decompo QR.
 */
