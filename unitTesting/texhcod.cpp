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




void printRand( int r,int c)
{

  std::cout << "<< ";
  for( int i=0;i<r;++i )
    {
      for( int j=0;j<c;++j )
	{
	  std::cout << std::setprecision(3) << ((rand()+0.0)/RAND_MAX*2)-1.;
	  if( j==c-1 )
	    { if( i==r-1 ) std::cout << " ;"; else std::cout << " ,"; }
	  else std::cout << " ,  ";
	}
      if( i==r-1 ) std::cout << std::endl << std::endl;
      else std::cout << std::endl << "         ";
    }
}

#define PRINT_RAND( A ) std::cout << "      " << #A ; printRand(A.rows(),A.cols());

void printStage( int s,int i,int r, int j )
{
  Eigen::MatrixXd Xhi( i,r );
  Eigen::MatrixXd Jfr( r,j );

  std::cout << "    {\n"
    "      s = "<<s<<"; // Stage s -- size "<<i<<","<<j<<" -- rank "<<r<<"\n"
    "      Eigen::MatrixXd Xhi( NR[ s],RANK[ s] );\n"
    "      Eigen::MatrixXd Jfr( RANK[ s],NC );\n"
    "      b[ s].resize(NR[ s]);\n" << std::endl;
  PRINT_RAND(Xhi); PRINT_RAND(Jfr);
  std::cout <<"      J[s] = Xhi*Jfr;\n"
    "      for( unsigned int i=0;i<NR[s];++i ) b[s][i] = (double)(i+1);\n"
    "    }"<<std::endl;
}

    // /* Automatic code generation.*/
    // if(0)
    //   {
    // 	for( int s=0;s<NB_STAGE;++s )
    // 	  printStage(s,NR[s],RANK[s],NC );
    // 	exit(0);
    //   }







void generateDataSet( std::vector<Eigen::MatrixXd> &J,
		      std::vector<soth::bound_vector_t> &b,
		      soth::HCOD& hcod,
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

  /* SOTH structure construction. */
  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      hcod.pushBackStage( J[i],b[i] );
      hcod.setInitialActiveSet( Eigen::VectorXi::LinSpaced(0, NR[i]-1, NR[i]),i);
    }
}


int main (int argc, char** argv)
{
  sotDebugTrace::openFile();
  bool exitOk=true;

  Eigen::MatrixXd Massert(5,9);
  soth::MatrixRnd::randomize( Massert );
  assert( std::abs(Massert(1,2)-0.985007)<1e-5 );

  {
    /* All matrices full rank, updating rows until the last rank saturates
     * du to the matrix size.
     */
    const int NB_STAGE = 3;
    const int RANK[] = { 3,3,6 };
    const int NR[] = { 3,3,6 };
    const int NC = 10;

    soth::HCOD hcod(NC,NB_STAGE);
    std::vector<Eigen::MatrixXd> J(NB_STAGE);
    std::vector<soth::bound_vector_t> b(NB_STAGE);
    generateDataSet( J,b,hcod, NB_STAGE,RANK,NR,NC );

    hcod.initialize();
    exitOk&=hcod.testRecomposition(&std::cout);
  }

  exit( (exitOk?0:1) );

}
