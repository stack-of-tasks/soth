/*
 *  Copyright
 */
#define SOTH_DEBUG
#define SOTH_DEBUG_MODE 45
#include "soth/debug.hpp"
#include "soth/HCOD.hpp"
#include "MatrixRnd.hpp"
#include <sys/time.h>
#include <fstream>
#include <Eigen/SVD>

using namespace soth;
using std::endl;
using std::cout;
using std::cerr;


void generateDeficientDataSet( std::vector<Eigen::MatrixXd> &J,
			       std::vector<soth::bound_vector_t> &b,
			       const int NB_STAGE,
			       const std::vector<int> & RANKFREE,
			       const std::vector<int> & RANKLINKED,
			       const std::vector<int> & NR,
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

      for( unsigned int i=0;i<NR[s];++i )
	{
	  double x = Random::rand<double>() * 2; // DEBUG: should b U*2-1
	  double y = Random::rand<double>() * -2;  // DEBUG
	  switch( randu(1,4) )
	    {
	    case 1: // =
	      b[s][i] = x;
	      break;
	    case 2: // <
	      b[s][i] = soth::Bound(y,soth::Bound::BOUND_INF);
	      break;
	    case 3: // >
	      b[s][i] = soth::Bound(x,soth::Bound::BOUND_SUP);
	      break;
	    case 4: // < <
	      b[s][i] = std::make_pair( std::min(x,y),std::max(x,y) ) ;
	      break;
	    }
	  if(s==3)
	    {
	      // if( i==0 ) b[s][i]=x;
	      // if( i==1 ) b[s][i]=y;
	      // if( i==2 ) b[s][i]=x;
	    }
	}
    }
}

void generateRandomProfile(int & nbStage,
			   std::vector<int>& rankfree,
			   std::vector<int>& ranklinked,
			   std::vector<int>& nr,
			   int & nc )
{
  nc = Random::rand<int>() % 50 + 6;
  nbStage = randu(1,1+nc/5);

  sotDEBUG(1) << "nc = " << nc << endl;
  sotDEBUG(1) << "nbStage = " << nbStage << endl;

  const int NR = std::max(2,(int)round((0.+nc)/nbStage*.7));
  const int RANKFREE = std::max(1,(int)round(whiteNoise(NR/2,0.2)));
  const int RANKLINKED = round(whiteNoise(NR,1))+1;
  sotDEBUG(1) << "mean_NR = " << NR << "; mean_RF = " << RANKFREE << "; mean_RL = " << RANKLINKED << endl;

  rankfree.resize( nbStage );
  ranklinked.resize( nbStage );
  nr.resize( nbStage );
  for( int i=0;i<nbStage;++i )
    {
      if( Random::rand<double>()<0.7 )
	{
	  sotDEBUG(1) << i<<": normal rank." <<endl;
	  rankfree[i] = randu(0,RANKFREE);//whiteNoise( RANKFREE,3 );
	  ranklinked[i] = randu(0,RANKLINKED); //whiteNoise( RANKLINKED,3 );
	  nr[i] = randu(1,NR);
	}
      else if( Random::rand<double>()<0.05 )
	{
	  sotDEBUG(1) << i<<":  rank def." <<endl;
	  rankfree[i] = randu(0,RANKFREE);//whiteNoise( RANKFREE,3 );
	  ranklinked[i] = randu(0,RANKLINKED); //whiteNoise( RANKLINKED,3 );
	  nr[i] = randu(1,NR);
	}
      else
	{
	  sotDEBUG(1) << i<<": full rank." <<endl;
	  ranklinked[i] = whiteNoise( RANKLINKED,3 );
	  nr[i] = randu(1,NR);
	  rankfree[i] = nr[i];
	}
      rankfree[i]=std::min(nr[i],rankfree[i]);
      ranklinked[i]=std::min(nr[i],ranklinked[i]);
      if( i==0 ) { ranklinked[i] = 0; rankfree[i]=std::max(1,rankfree[i] ); }
      else ranklinked[i]=std::min( ranklinked[i],nr[i-1] );
      if( rankfree[i]==0 ) ranklinked[i]=std::max(1,ranklinked[i]);
      sotDEBUG(1) << "rf"<<i<<" = " << rankfree[i] <<";   rl"<<i<<" = " << ranklinked[i]
		  << ";  nr"<<i<<" = " << nr[i] << endl;

    }
}


/* Compute [m1;m2] from m1 and m2 of same col number. */
template< typename D1,typename D2 >
MatrixXd stack( const MatrixBase<D1>& m1, const MatrixBase<D2>& m2 )
{
  assert( m1.cols() == m2.cols() );
  const int m1r=m1.rows(), m2r=m2.rows(), mr=m1r+m2r, mc=m1.cols();
  MatrixXd res( mr,mc );
  for( int i=0;i<m1r;++i )
    for( int j=0;j<mc;++j ) res(i,j) = m1(i,j);
  for( int i=0;i<m2r;++i )
    for( int j=0;j<mc;++j ) res(m1r+i,j) = m2(i,j);
  return res;
}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
int main (int argc, char** argv)
{
  {
    struct timeval tv;
    gettimeofday(&tv,NULL);

    int seed = tv.tv_usec % 7919; //= 7594;
    if( argc == 2 )
      {  seed = atoi(argv[1]);  }
    std::cout << "seed = " << seed << std::endl;
    soth::Random::setSeed(seed);
  }

  /* Decide the size of the problem. */
  int NB_STAGE,NC;
  std::vector<int> NR,RANKLINKED,RANKFREE;
  generateRandomProfile(NB_STAGE,RANKFREE,RANKLINKED,NR,NC);

  /* Initialize J and b. */
  std::vector<Eigen::MatrixXd> J(NB_STAGE);
  std::vector<soth::bound_vector_t> b(NB_STAGE);
  generateDeficientDataSet(J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC);

  std::ofstream fout("/tmp/soth.txt");
  fout << "variable size " << NC << endl << endl;

  for( unsigned int s=0;s<NB_STAGE;++s )
    {
      fout << "level" << endl << endl;

      int nbEqualities = 0;
      for( unsigned int r=0;r<NR[s];++r )
	{
	  if( b[s][r].getType() == Bound::BOUND_TWIN ) nbEqualities++;
	}

      fout << "equalities " << nbEqualities << endl << endl;
      for( unsigned int r=0;r<NR[s];++r )
	{
	  if( b[s][r].getType() == Bound::BOUND_TWIN )
	    {
	      fout << J[s].row(r) << "   " << b[s][r].getBound( Bound::BOUND_TWIN ) << endl;
	    }
	}
      fout << endl << "inequalities " << NR[s]-nbEqualities << endl << endl;
      for( unsigned int r=0;r<NR[s];++r )
	{
	  switch(b[s][r].getType())
	    {
	    case Bound::BOUND_TWIN:
	      break;
	    case Bound::BOUND_INF:
	      fout << J[s].row(r) << "   " << b[s][r].getBound( Bound::BOUND_INF ) << " 1e25" << endl;
	      break;
	    case Bound::BOUND_SUP:
	      fout << J[s].row(r) << "   " << -1e25 << " " << b[s][r].getBound( Bound::BOUND_SUP ) << endl;
	      break;
	    case Bound::BOUND_DOUBLE:
	      fout << J[s].row(r) << "   " <<  b[s][r].getBound( Bound::BOUND_INF ) << " " << b[s][r].getBound( Bound::BOUND_SUP ) << endl;
	      break;
	    }
	}
      fout << endl;
    }

  fout << endl << "end" << endl << endl;
}
