/*
 *  Copyright
 */
#define SOTH_DEBUG
#define SOTH_DEBUG_MODE 45
#include "soth/debug.h"
#include "soth/HCOD.hpp"
#include "soth/debug.h"
#include "MatrixRnd.h"
#include <sys/time.h>
#include <iostream>
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

/* --------------------------------------------------------------------------- */
/* --- COD SOLVER ------------------------------------------------------------ */
/* --------------------------------------------------------------------------- */
/* Compute the decomposition A=U.L.V', with U and L rotation matrices, and
 * L = [ L0 0 ; 0 0 ] with L0 triangular inf, with non zero diagonal.
 * Use this decompo to solve min||Ax-b||.
 *
 * This class is for debug only. The computations involved in the solver
 * along with the ones in the projection are very suboptimal, and should
 * be properly rewritten for real times used.
 */
struct ULV
{
  ColPivHouseholderQR<MatrixXd> qru;
  MatrixXd Rt,Ainit;
  ColPivHouseholderQR<MatrixXd> qrv;
  int rank,NC,NR;

  void compute( const MatrixXd& A, const unsigned int & rank_ )
  {
    rank=rank_; NC=A.cols(); NR=A.rows();
#ifdef DEBUG
    Ainit=A;
#endif

    if( rank==0 ) return;

    qru.compute( A );
    Rt = qru.matrixQR().topRows(rank).triangularView<Upper>().transpose();
    qrv.compute( Rt );
  }

  void disp( bool check=false )
  {
    if( rank==0 )
      {
	sotDEBUG(5) << "ULV: Empty rank."  << endl;
	return;
      }

    sotDEBUG(5) << "U = " << (MATLAB)(MatrixXd)qru.householderQ() << endl;
    sotDEBUG(5) << "R = " << (MATLAB)(MatrixXd)Rt.transpose() << endl;
    sotDEBUG(5) << "V = " << (MATLAB)(MatrixXd)qrv.householderQ() << endl;
    sotDEBUG(5) << "L = " << (MATLAB)(MatrixXd)qrv.matrixQR().topRows(rank).triangularView<Upper>().transpose() << endl;
    sotDEBUG(5) << "Pcv = " << (MATLAB)(MatrixXd)(qru.colsPermutation()) << endl;
    sotDEBUG(5) << "Pcu = " << (MATLAB)(MatrixXd)(qrv.colsPermutation()) << endl;

    if( check )
      {
	MatrixXd L0(NR,NC); L0.setZero();
	L0.topLeftCorner(rank,rank)
	  = qrv.matrixQR().topRows(rank).triangularView<Upper>().transpose();
	MatrixXd Arec = matrixU() * L0 * matrixV().transpose();
	sotDEBUG(5) << "Arec = " << (MATLAB)Arec << endl;
	assert( (Ainit-Arec).norm()< Stage::EPSILON );
      }
  }

  /* Solve min||Ax-b|| for a matrix A whose rank is given. */
  VectorXd solve( const VectorXd& b ) const
  {
    if( rank==0 ) return VectorXd::Zero(NC);

    VectorXd sol = b;  /* s = b */
    sol.applyOnTheLeft(qru.householderQ().adjoint()); /* s = U'*s */
    sol.applyOnTheLeft(qrv.colsPermutation().transpose());  /* s = Pu'*s */
    sotDEBUG(5) << "PuUtb = " << (MATLAB)sol;

    VectorXd solv(NC); solv.setZero();
    solv.head(rank) = sol.head(rank);
    qrv.matrixQR().topRows(rank).transpose().triangularView<Lower>()
      .solveInPlace( solv.head(rank) );  /* s = Linv*s */
    sotDEBUG(5) << "LiPuUtb = " << (MATLAB)solv;

    solv.applyOnTheLeft(qrv.householderQ());  /* s = V*s */
    solv.applyOnTheLeft(qru.colsPermutation()); /* s = Pv*s */
    sotDEBUG(5) << "VLUe = " << (MATLAB)solv << endl;
    return solv;
  }

  MatrixXd matrixV(void) const
  {
    if( rank==0 ){ return MatrixXd::Identity(NC,NC); }
    MatrixXd V = qru.colsPermutation();
    V.applyOnTheRight(qrv.householderQ());
    return V;
  }
  MatrixXd matrixU(void) const
  {
    if( rank==0 ){ return MatrixXd::Identity(NR,NR); }
    MatrixXd U = MatrixXd::Identity(NR,NR);
    U.topLeftCorner(rank,rank) = qrv.colsPermutation(); /* I think Pu is always 1.*/
    U.applyOnTheLeft(qru.householderQ());  /* U = H*U */
    return U;
  }
  void decreaseProjector( MatrixXd & P ) const
  { /* Highly suboptimal ... */
    if( rank==0 ) return;
    MatrixXd V1 = matrixV().leftCols(rank);
    P -= V1*V1.transpose();
  }
};

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
  sotDebugTrace::openFile();

  /* Decide the size of the problem. */
  int NB_STAGE,NC;
  std::vector<int> NR,RANKLINKED,RANKFREE;
  generateRandomProfile(NB_STAGE,RANKFREE,RANKLINKED,NR,NC);

  /* Initialize J and b. */
  std::vector<Eigen::MatrixXd> J(NB_STAGE);
  std::vector<soth::bound_vector_t> b(NB_STAGE);
  generateDeficientDataSet(J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC);

  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      std::cout << "J"<<i+1<<" = " << (soth::MATLAB)J[i] << std::endl;
      std::cout << "e"<<i+1<< " = " << b[i] << ";"<<std::endl;
    }
  //  assert( std::abs(J[0](0,0)-(-1.1149))<1e-5 );

  /* SOTH structure construction. */
  soth::HCOD hcod(NC,NB_STAGE);
  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      hcod.pushBackStage( J[i],b[i] );
    }
  hcod.setNameByOrder("stage_");

  VectorXd solution;
  hcod.activeSearch( solution );
  if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );


  /* --- CHECK --- */
  VectorXd u=solution,du = VectorXd::Zero(NC);
  MatrixXd Ja,Japrec;
  VectorXd ea,eaprec;
  MatrixXd Pa = MatrixXd::Identity(NC,NC);
  VectorXd usvd = VectorXd::Zero(NC);
  const double EPSILON = 10*Stage::EPSILON;

  for( unsigned int i=0;i<hcod.nbStages();++i )
    {
      Stage & st = hcod[i];
      MatrixXd J_(NR[i],NC); VectorXd e_(NR[i]);

      std::cout << "Stage " << i << "... " << endl;
      {
	for( int r=0;r<NR[i];++r )
	  {
	    double Ju = J[i].row(r)*solution;
	    cout << "   "<<i << ":" << r << "\t";
	    switch( b[i][r].getType() )
	      {
	      case Bound::BOUND_TWIN:
		{
		  double x = b[i][r].getBound( Bound::BOUND_TWIN );
		  if( std::abs( x-Ju )<EPSILON )
		    cout << " (=)    :  \t  Ju="  << Ju
			 << " ~ " << x << "=b" << endl;
		  else
		    cerr << "!! " << " (=):  \t  Ju="  << Ju
			 << " != " << x << "=b" << endl;
		  break;
		}
	      case Bound::BOUND_INF:
		{
		  double x = b[i][r].getBound( Bound::BOUND_INF );
		  if( x<Ju+EPSILON )
		    cout << " (b<.)  :  \t  Ju="  << Ju
			 << " > " << x << "=b" << endl;
		  else
		    cerr << "!! " << " (=):  \t  Ju="  << Ju
			 << " !< " << x << "=b" << endl;
		  break;
		}
	      case Bound::BOUND_SUP:
		{
		  double x = b[i][r].getBound( Bound::BOUND_SUP );
		  if( Ju<x+EPSILON )
		    cout << " (.<b)  :  \t  Ju="  << Ju
			 << " < " << x << "=b" << endl;
		  else
		    cerr << "!! " << " (=):  \t  Ju="  << Ju
			 << " !> " << x << "=b" << endl;
		  break;
		}
	      case Bound::BOUND_DOUBLE:
		{
		  double xi = b[i][r].getBound( Bound::BOUND_INF );
		  double xs = b[i][r].getBound( Bound::BOUND_SUP );
		  if( (xi<Ju+EPSILON)&&(Ju<xs+EPSILON) )
		    cout << " (b<.<b):  \t  binf="<<xi<<" < Ju="  << Ju
			 << " < " << xs << "=bsup" << endl;
		  else
		    cout << "!! (b<.<b):  \t  binf="<<xi<<" !< Ju="  << Ju
			 << " !< " << xs << "=bsup" << endl;
		  break;
		}
	      }


	  }
      }


      sotDEBUG(1) << "Check bounds of " << i << "."<<endl;
      assert( st.checkBound(u,du,NULL,NULL) );

      if( st.sizeA()==0 ) continue;

      SubMatrix<MatrixXd,RowPermutation> Jai = st.Jactive(J_);
      SubMatrix<VectorXd,RowPermutation> eai = st.eactive(e_);

      if( st.rank() == st.sizeA() )
	{
	  sotDEBUG(1) << "Check fullrankness of " << i << "."<<endl;
	  assert( ( eai - Jai*u ).norm()<st.EPSILON );
	}

      MatrixXd JPi = Jai*Pa;
      const int rank = st.rank();
      sotDEBUG(5) << "e"<<i<< " = " << (MATLAB)eai << endl;
      sotDEBUG(5) << "J"<<i<< " = " << (MATLAB)Jai << endl;
      sotDEBUG(5) << "JP"<<i<< " = " << (MATLAB)JPi << endl;
      sotDEBUG(5) << "rank"<<i<< " = " << (MATLAB)rank << endl;

      ULV ulv; ulv.compute(JPi,rank); usvd += ulv.solve( eai-Jai*usvd );
      ulv.decreaseProjector( Pa );

      ulv.disp(true);
      sotDEBUG(5) << "usvd"<<i<< " = " << (MATLAB)usvd << endl;
      sotDEBUG(5) << "V"<<i<< " = " << (MATLAB)ulv.matrixV() << endl;
      sotDEBUG(5) << "P"<<i<< " = " << (MATLAB)Pa << endl;

      sotDEBUG(1) << "Check pinv of " << i << "." << endl;
      assert( std::abs(( eai-Jai*u ).norm() - ( eai-Jai*usvd ).norm()) < 10*Stage::EPSILON );
    }

}