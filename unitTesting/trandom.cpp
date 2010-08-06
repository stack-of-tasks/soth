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

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#include <fstream>
void parseProblemFile( const std::string name,
		       std::vector<Eigen::MatrixXd> &J,
		       std::vector<soth::bound_vector_t> &b,
		       int& NB_STAGE,
		       std::vector<int> & NR,
		       int& NC )
{
  std::ifstream fin(name.c_str());
  std::string syntax1,syntax2;
  MatrixXd Ji,eiinf,eisup;

  fin >> syntax1 >> syntax2 >> NC;
  assert( (syntax1 == "variable")&&(syntax2 == "size") );
  NB_STAGE=0;  int & s= NB_STAGE;
  fin >> syntax1;
  do
    {
      NR.resize(s+1); J.resize(s+1); b.resize(s+1);

      unsigned int nre;
      fin >> syntax2 >> nre;
      assert( (syntax1 == "level")&&(syntax2 == "equalities") );

      MatrixXd Je; VectorXd ee;
      if( nre>0 )
	{
	  Je.resize(nre,NC); ee.resize(nre);
	  for( unsigned int i=0;i<nre;++i )
	    {
	      for( unsigned int j=0;j<NC;++j )
		fin >> Je(i,j);
	      fin >> ee(i);
	    }
	}

      fin >> syntax1 >>  NR[s];
      assert( (syntax1 == "inequalities") );

      /* Copy the equalities line into the total matrix. */
      NR[s]+=nre; assert(NR[s]>0);
      J[s].resize(NR[s],NC); b[s].resize(NR[s]);
      if( nre>0 )
	{
	  J[s].topRows(nre)=Je;
	  for( unsigned int i=0;i<nre;++i )  b[s][i] = ee(i);
	}

      /* Parse the inequalities. */
      if( NR[s]>nre )
	{
	  double bi,bu;
	  for( unsigned int i=nre;i<NR[s];++i )
	    {
	      for( unsigned int j=0;j<NC;++j )
		fin >> J[s](i,j);
	      fin >> bi>>bu;
	      if( bi<-1e10 ) // bound sup only
		{
		  assert( bu<=1e10 );
		  b[s][i] = Bound( bu,Bound::BOUND_SUP );
		}
	      else if( 1e10<bu ) // bound inf only
		{
		  b[s][i] = Bound( bi,Bound::BOUND_INF );
		}
	      else // double bound
		{
		  b[s][i] = std::make_pair( bi,bu );
		}
	    }
	}
      sotDEBUG(5) << "J"<<s<<" = " << (MATLAB)J[s] << endl;
      sotDEBUG(5) << "b"<<s<<" = " << b[s] << endl;
      fin >> syntax1; s++;
    } while( syntax1 == "level" );

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
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
	  int btype = randu(1,4);
	  if( i==0 ) btype= randu(2,4);
	  switch( btype )
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
  nbStage = 2; //randu(1,1+nc/5);

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

  void compute( const MatrixXd& A, const double & svmin )
  {
    NC=A.cols(); NR=A.rows();
    const int N=std::min(NC,NR);
#ifdef DEBUG
    Ainit=A;
#endif

    qru.compute( A );

    const MatrixXd& QR = qru.matrixQR();
    for( rank=0;rank<N;++rank ) if( std::abs(QR(rank,rank))<=svmin ) break;
    if( rank==0 ) return;

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
/* --- SOTH SOLVER ---------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* This piece of code tries all the possible active set, looking for the best
 * solution, and compare this solution with the current solution. The function
 * to be run is explore, with the constraints J and b as arguments, as well
 * as the optimum to be checked. */
namespace DummyActiveSet
{


  VectorXd SOT_Solver( std::vector<Eigen::MatrixXd> J,
		       std::vector<VectorXd> e )
  {
    const unsigned int NC = J[0].cols();

    VectorXd u = VectorXd::Zero(NC);
    MatrixXd Pa = MatrixXd::Identity(NC,NC);
    const double EPSILON = Stage::EPSILON;

    for( unsigned int i=0;i<J.size();++i )
      {
	if( J[i].rows() == 0 ) continue;
	MatrixXd JPi = J[i]*Pa;
	ULV ulv; ulv.compute(JPi,EPSILON);
	u += ulv.solve( e[i]-J[i]*u );
	ulv.decreaseProjector( Pa );
      }

    return u;


  }

  std::vector<double> stageErrors( const std::vector<Eigen::MatrixXd>& Jref,
				   const std::vector<soth::bound_vector_t>& bref,
				   const VectorXd& usot )
  {
    std::vector<double> errors(Jref.size());
    for( unsigned int i=0;i<Jref.size();++i )
      {
	const MatrixXd& J = Jref[i];
	const soth::bound_vector_t& b = bref[i];

	VectorXd e = J*usot;
	errors[i]=0;
	for( unsigned int r = 0;r<J.rows();++r )
	  {
	    double x = b[r].distance(e(r));
	    errors[i] += x*x;
	  }
	sotDEBUG(5) << "err(" << i << ") = " << errors[i] << endl;
      }
    return errors;
  }

  /* Alphalexical order */
  //template<>
  bool isLess ( const std::vector< double > & m1, const std::vector<double> & m2 )
  {
    assert( m1.size()==m2.size() );
    for( unsigned int i=0;i<m1.size();++i )
      {
	if( m1[i]-Stage::EPSILON>m2[i] ) return false;
      }
    return true;
  }


  /* Use the old SOT solver for the given active set. The optimum
   * is then compared to the bound. If the found solution is valid, compare its
   * norm to the reference optimum. Return true is the solution found by SOT is valid and
   * better. */
  bool compareSolvers( const std::vector<Eigen::MatrixXd>& Jref,
		       const std::vector<soth::bound_vector_t>& bref,
		       const std::vector< std::vector<int> >& active,
		       const std::vector< std::vector<Bound::bound_t> >& bounds,
		       const double & urefnorm,
		       const std::vector<double> & erefnorm,
		       const bool verbose = false )
  {
    /* Build the SOT problem. */
    const unsigned int NC = Jref[0].cols();
    std::vector<Eigen::MatrixXd> Jsot(Jref.size());
    std::vector<VectorXd> esot(Jref.size());
    for( unsigned int i=0;i<Jref.size();++i )
      {
	MatrixXd& J = Jsot[i]; VectorXd& e= esot[i];
	const std::vector<int> & Aset = active[i];
	const std::vector<Bound::bound_t> & Bset = bounds[i];
	const unsigned int NR = Aset.size();

	J.resize(NR,NC); e.resize(NR);
	for( unsigned int r=0;r<NR;++r )
	  {
	    sotDEBUG(5) << "Get " << i<<"," << Aset[r] <<" ("<<r<<")" << endl;
	    J.row(r) = Jref[i].row( Aset[r] );
	    e(r) = bref[i][Aset[r]].getBound( Bset[Aset[r]] );
	  }

	sotDEBUG(45) << "J" << i << " = " << (MATLAB)J << endl;
	sotDEBUG(45) << "e" << i << " = " << (MATLAB)e << endl;
      }

    /* Find the optimum. */
    VectorXd usot = SOT_Solver( Jsot,esot );
    sotDEBUG(45) << "usot" <<" = " << (MATLAB)usot << endl;
    //if( (!verbose)&&(usot.norm()+Stage::EPSILON >= urefnorm) ) return false;

    /* Check the bounds. */
    std::vector<double> esotnorm = stageErrors(Jref,bref,usot);
    if( verbose||isLess( esotnorm,erefnorm ) )
      {
	std::cout << endl << endl << " * * * * * * * * * * * * * " << endl
		  <<  "usot" <<" = " << (MATLAB)usot << endl;
	std::cout << "aset = {";
	for( unsigned int i=0;i<Jref.size();++i )
	  {
	    std::cout << "\n        [ ";
	    for( unsigned int r=0;r<active[i].size();++r )
	      {
		const int ref = active[i][r];
		if( bref[i][ref].getType() ==  Bound::BOUND_DOUBLE )
		  if( bounds[i][ref] == Bound::BOUND_INF ) cout << "-";
		  else if( bounds[i][ref] == Bound::BOUND_SUP ) cout << "+";
		std::cout <<active[i][r] << " ";
	      }
	    std::cout << " ]";
	  }
	std::cout << "} " << endl;
	std::cout << "norm solution = " << usot.norm() << " ~?~ " << urefnorm << " = norm ref." << endl;
	std::cout << "e = [" ;
	for( unsigned int s=0;s<esotnorm.size();++s )
	  { std::cout << esotnorm[s] << "    "; }
	std::cout << " ];" << endl;

	return true;
      }
    else return false;

  }

  void intToVbool( const int nbConstraint, const unsigned long int ref,
		   std::vector<bool>& res )
  {
    res.resize(nbConstraint);
    sotDEBUG(45) << "ref" <<" = " << ref << endl;
    for( unsigned int i=0;i<nbConstraint;++i )
      {
	sotDEBUG(45) << i << ": " << (0x01&(ref>>i)) << endl;
	res[i]=( 0x01&(ref>>i) );
      }
  }

  void selectConstraint( const std::vector<soth::bound_vector_t>& bref,
			 const std::vector<bool>& activeBool,
			 std::vector< std::vector<int> > & aset )
  {
    aset.resize(bref.size());
    int cst = 0;
    for( unsigned int i=0;i<bref.size();++i )
      {
	for( unsigned int r=0;r<bref[i].size();++r )
	  if( (bref[i][r].getType() == Bound::BOUND_TWIN)
	      || (activeBool[cst++]) ) aset[i].push_back(r);

	// std::cout << "a" << i << " = [ ";
	// for( unsigned int r=0;r<aset[i].size();++r )
	// 	std::cout << aset[i][r] << " ";
	// std::cout << "];" << endl;
      }
  }

  void selectBool( const std::vector<soth::bound_vector_t>& bref,
		   std::vector<std::vector<int> > aset,
		   const std::vector<bool>& boundBool,
		   std::vector< std::vector<Bound::bound_t> >& boundSelec )
  {
    int boolPrec=0;
    boundSelec.resize(bref.size());
    for( unsigned int i=0;i<bref.size();++i )
      {
	boundSelec[i].resize(bref[i].size(),Bound::BOUND_NONE);
	for( unsigned int r=0;r<aset[i].size();++r )
	  {
	    if( bref[i][aset[i][r]].getType() == Bound::BOUND_DOUBLE )
	      {
		sotDEBUG(5) << "Fill " << i<<"," << aset[i][r] <<" ("<<r<<")" << endl;
		boundSelec[i][aset[i][r]] = (boundBool[boolPrec++])?Bound::BOUND_INF:Bound::BOUND_SUP;
	      }
	    else
	      {
		boundSelec[i][aset[i][r]] = bref[i][aset[i][r]].getType();
	      }
	    assert( (boundSelec[i][aset[i][r]]!=Bound::BOUND_DOUBLE)&&(boundSelec[i][aset[i][r]]!=Bound::BOUND_NONE) );
	  }
      }
  }

  int computeNbDouble( const std::vector<soth::bound_vector_t>& bref,
		       std::vector<std::vector<int> > aset )
  {
    int nbDoubleConstraint = 0;
    for( unsigned int i=0;i<bref.size();++i )
      {
	const soth::bound_vector_t& b = bref[i];
	for( unsigned int r = 0;r<aset[i].size();++r )
	  {
	    if( bref[i][aset[i][r]].getType() == Bound::BOUND_DOUBLE ) nbDoubleConstraint++;
	  }
      }
    return nbDoubleConstraint;
  }

  int computeNbConst( const std::vector<soth::bound_vector_t>& bref,
		      std::vector<int>& nr )
  {
    int nbConstraint = 0; nr.resize(bref.size());
    for( unsigned int i=0;i<bref.size();++i )
      {
	nr[i]=bref[i].size();
	for( unsigned int r = 0;r<bref[i].size();++r )
	  if(bref[i][r].getType() != Bound::BOUND_TWIN)
	    nbConstraint ++;
      }
    return nbConstraint;
  }


  /* Explore all the possible active set.
   */
  void explore( const std::vector<Eigen::MatrixXd>& Jref,
		const std::vector<soth::bound_vector_t>& bref,
		const VectorXd& uref )
  {
    std::vector<int> nr(Jref.size());
    int nbConstraint = computeNbConst(bref,nr);
    const double urefnorm = uref.norm();
    const std::vector<double> erefnorm = stageErrors( Jref,bref,uref );

    std::cout << "Nb posibilities = " << (1<<(nbConstraint)) << endl;
    for( unsigned long int refc =0;refc< (1<<(nbConstraint)); ++refc )
      {
	if( !( refc%1000) ) std::cout << refc << " ... " << endl;
	std::vector<bool> abool; std::vector<std::vector<int> > aset;
	intToVbool(nbConstraint,refc,abool); selectConstraint(bref,abool,aset);

	int nbDoubleConstraint = computeNbDouble(bref,aset);
	sotDEBUG(15) << "Nb double = " << nbDoubleConstraint << endl;
	for( unsigned long int refb =0;refb< (1<<(nbDoubleConstraint)); ++refb )
	  {
	    std::vector<bool> bbool; std::vector< std::vector<Bound::bound_t> > bset;

	    intToVbool(nbDoubleConstraint,refb,bbool); selectBool( bref,aset,bbool,bset );

	    if(  compareSolvers( Jref,bref,aset,bset,urefnorm,erefnorm ) )
	      { std::cout << refc << "/" << refb << "  ...  One more-optimal solution!! " << endl; }
	  }
      }
  }


  void detailActiveSet( const std::vector<Eigen::MatrixXd>& Jref,
			const std::vector<soth::bound_vector_t>& bref,
			const VectorXd& uref,
			unsigned long int refc,unsigned long int refb )
  {
    std::vector<int> nr(Jref.size());
    int nbConstraint = computeNbConst(bref,nr);
    const double urefnorm = uref.norm();
    const std::vector<double> erefnorm = stageErrors( Jref,bref,uref );

    std::vector<bool> abool; std::vector<std::vector<int> > aset;
    intToVbool(nbConstraint,refc,abool); selectConstraint(bref,abool,aset);
    int nbDoubleConstraint = computeNbDouble(bref,aset);
    std::vector<bool> bbool; std::vector< std::vector<Bound::bound_t> > bset;
    intToVbool(nbDoubleConstraint,refb,bbool); selectBool( bref,aset,bbool,bset );

    compareSolvers( Jref,bref,aset,bset,urefnorm,erefnorm,true );
  }
};



/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
int main (int argc, char** argv)
{
  sotDebugTrace::openFile();

  int NB_STAGE,NC;
  std::vector<int> NR,RANKLINKED,RANKFREE;
  std::vector<Eigen::MatrixXd> J;
  std::vector<soth::bound_vector_t> b;

  if( (argc==3)&& std::string(argv[1])=="-file")
    {
      parseProblemFile( argv[2],J,b,NB_STAGE,NR,NC);
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
  cout << "Optimal active set = "; hcod.showActiveSet(std::cout);

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
	    cout << "   "<<i << ":" << r << (hcod[i].isActive(r)?"a":" ") <<"\t";
	    switch( b[i][r].getType() )
	      {
	      case Bound::BOUND_TWIN:
		{
		  double x = b[i][r].getBound( Bound::BOUND_TWIN );
		  if( std::abs( x-Ju )<EPSILON )
		    cout << "    (=)    :  \t  Ju="  << Ju
			 << " ~ " << x << "=b" << endl;
		  else
		    cout << "!! " << " (=):  \t  Ju="  << Ju
			 << " != " << x << "=b" << endl;
		  break;
		}
	      case Bound::BOUND_INF:
		{
		  double x = b[i][r].getBound( Bound::BOUND_INF );
		  if( x<Ju+EPSILON )
		    cout << "    (b<.)  :  \t  Ju="  << Ju
			 << " > " << x << "=b" << endl;
		  else
		    cout << "!! " << " (b<.):  \t  Ju="  << Ju
			 << " !< " << x << "=b" << endl;
		  break;
		}
	      case Bound::BOUND_SUP:
		{
		  double x = b[i][r].getBound( Bound::BOUND_SUP );
		  if( Ju<x+EPSILON )
		    cout << "    (.<b)  :  \t  Ju="  << Ju
			 << " < " << x << "=b" << endl;
		  else
		    cout << "!! " << " (.<b):  \t  Ju="  << Ju
			 << " !> " << x << "=b" << endl;
		  break;
		}
	      case Bound::BOUND_DOUBLE:
		{
		  double xi = b[i][r].getBound( Bound::BOUND_INF );
		  double xs = b[i][r].getBound( Bound::BOUND_SUP );
		  if( (xi<Ju+EPSILON)&&(Ju<xs+EPSILON) )
		    cout << "    (b<.<b):  \t  binf="<<xi<<" < Ju="  << Ju
			 << " < " << xs << "=bsup" << endl;
		  else
		    cout << "!!  (b<.<b):  \t  binf="<<xi<<" !< Ju="  << Ju
			 << " !< " << xs << "=bsup" << endl;
		  break;
		}
	      }


	  }
      }


      sotDEBUG(1) << "Check bounds of " << i << "."<<endl;
      // DEBUG assert( st.checkBound(u,du,NULL,NULL) );

      if( st.sizeA()==0 ) continue;

      SubMatrix<MatrixXd,RowPermutation> Jai = st.Jactive(J_);
      SubMatrix<VectorXd,RowPermutation> eai = st.eactive(e_);

      if( st.rank() == st.sizeA() )
	{
	  sotDEBUG(1) << "Check fullrankness of " << i << "."<<endl;
	  assert( ( eai - Jai*u ).norm()<st.EPSILON );
	}

      MatrixXd JPi = Jai*Pa;
      const unsigned int rank = st.rank();
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



  DummyActiveSet::explore(J,b,solution);
  // DummyActiveSet::detailActiveSet( J,b,solution,2720,0 );
  // DummyActiveSet::detailActiveSet( J,b,solution,2721,1 );
  // DummyActiveSet::detailActiveSet( J,b,solution,2736,0 );
}
