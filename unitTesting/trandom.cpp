/*
 *  Copyright
 */
#define SOTH_DEBUG
#define SOTH_DEBUG_MODE 45
#include "soth/debug.hpp"
#include "soth/HCOD.hpp"
#include "soth/COD.hpp"
#include "MatrixRnd.hpp"
#include "RandomGenerator.hpp"
#include <sys/time.h>
#include <iostream>

using namespace soth;
using std::endl;
using std::cout;
using std::cerr;

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */


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
	for( unsigned int r = 0;int(r)<J.rows();++r )
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
	if( m1[i]>m2[i] ) return false;
      }
    return true;
  }

  std::vector< double >
  operator+ ( const std::vector< double > & m1, const double & b)
  {
    std::vector< double > res; res.resize(m1.size() );
    for( unsigned int i=0;i<m1.size();++i ) res[i]=m1[i]+b;
    return res;
  }
  std::vector< double >
  operator- ( const std::vector< double > & m1, const double & b)
  {
    std::vector< double > res; res.resize(m1.size() );
    for( unsigned int i=0;i<m1.size();++i ) res[i]=m1[i]-b;
    return res;
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
    if( verbose
	||(isLess( esotnorm+Stage::EPSILON,erefnorm ) ) )
	  //	   &&( usot.norm()+Stage::EPSILON < urefnorm) ) ) //DEBUG
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
		if( (bref[i][ref].getType()) ==  Bound::BOUND_DOUBLE )
		  {
		    if( bounds[i][ref] == Bound::BOUND_INF ) cout << "-";
		    else if( bounds[i][ref] == Bound::BOUND_SUP ) cout << "+";
		  }
		std::cout <<active[i][r] << " ";
	      }
	    std::cout << " ]";
	  }
	std::cout << "} " << endl;
	std::cout << "norm solution = " << usot.norm() << " ~?~ " << urefnorm << " = norm ref." << endl;
	std::cout << "e = [" ;
	for( unsigned int s=0;s<esotnorm.size();++s )
	  { std::cout << (esotnorm+Stage::EPSILON)[s] << "    (" << erefnorm[s] << ")   "; }
	std::cout << " ];" << endl;

	return true;
      }
    else return false;

  }

  void intToVbool( const unsigned int nbConstraint, const unsigned long int ref,
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
  bool explore( const std::vector<Eigen::MatrixXd>& Jref,
		const std::vector<soth::bound_vector_t>& bref,
		const VectorXd& uref )
  {
    std::vector<int> nr(Jref.size());
    int nbConstraint = computeNbConst(bref,nr);
    const double urefnorm = uref.norm();
    const std::vector<double> erefnorm = stageErrors( Jref,bref,uref );

    bool res=true;
    std::cout << "Nb posibilities = " << (1<<(nbConstraint)) << endl;
    for( long int refc =0;refc< (1<<(nbConstraint)); ++refc )
      {
	if( !( refc%1000) ) std::cout << refc << " ... " << endl;
	std::vector<bool> abool; std::vector<std::vector<int> > aset;
	intToVbool(nbConstraint,refc,abool); selectConstraint(bref,abool,aset);

	int nbDoubleConstraint = computeNbDouble(bref,aset);
	sotDEBUG(15) << "Nb double = " << nbDoubleConstraint << endl;
	for( long int refb =0;refb< (1<<(nbDoubleConstraint)); ++refb )
	  {
	    std::vector<bool> bbool; std::vector< std::vector<Bound::bound_t> > bset;

	    intToVbool(nbDoubleConstraint,refb,bbool); selectBool( bref,aset,bbool,bset );

	    if(  compareSolvers( Jref,bref,aset,bset,urefnorm,erefnorm ) )
	      {
		std::cout << refc << "/" << refb
			  << "  ...  One more-optimal solution!! " << endl;
		res=false;
	      }
	  }
      }
    return res;
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
# ifndef NDEBUG
#endif
  sotDebugTrace::openFile();

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
  cout << endl;
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

  //exit(0); // DEBUG

  /* --- CHECK --- */
  VectorXd u=solution,du = VectorXd::Zero(NC);
  MatrixXd Ja,Japrec;
  VectorXd ea,eaprec;
  MatrixXd Pa = MatrixXd::Identity(NC,NC);
  VectorXd usvd = VectorXd::Zero(NC);

  /* Checks and recheck ...
   * These checks are not valid: first, checking if the bound is
   * or is not respected has little meaning, when the stages are not full rank.
   * Second, when a slack is freezed, it is not possible to check its value any
   * more with the initial reference value. .. TODO again

     const double EPSILON = 10*Stage::EPSILON;
     for( unsigned int i=0;i<hcod.nbStages();++i )
    {
      Stage & st = hcod[i];
      MatrixXd J_(NR[i],NC); VectorXd e_(NR[i]);

      std::cout << "Stage " << i << "... " << endl;
      // {
      // 	for( int r=0;r<NR[i];++r )
      // 	  {
      // 	    double Ju = J[i].row(r)*solution;
      // 	    cout << "   "<<i << ":" << r << (hcod[i].isActive(r)?"a":" ") <<"\t";
      // 	    switch( b[i][r].getType() )
      // 	      {
      // 	      case Bound::BOUND_TWIN:
      // 		{
      // 		  double x = b[i][r].getBound( Bound::BOUND_TWIN );
      // 		  if( std::abs( x-Ju )<EPSILON )
      // 		    cout << "    (=)    :  \t  Ju="  << Ju
      // 			 << " ~ " << x << "=b" << endl;
      // 		  else
      // 		    cout << "!! " << " (=):  \t  Ju="  << Ju
      // 			 << " != " << x << "=b" << endl;
      // 		  break;
      // 		}
      // 	      case Bound::BOUND_INF:
      // 		{
      // 		  double x = b[i][r].getBound( Bound::BOUND_INF );
      // 		  if( x<Ju+EPSILON )
      // 		    cout << "    (b<.)  :  \t  Ju="  << Ju
      // 			 << " > " << x << "=b" << endl;
      // 		  else
      // 		    cout << "!! " << " (b<.):  \t  Ju="  << Ju
      // 			 << " !< " << x << "=b" << endl;
      // 		  break;
      // 		}
      // 	      case Bound::BOUND_SUP:
      // 		{
      // 		  double x = b[i][r].getBound( Bound::BOUND_SUP );
      // 		  if( Ju<x+EPSILON )
      // 		    cout << "    (.<b)  :  \t  Ju="  << Ju
      // 			 << " < " << x << "=b" << endl;
      // 		  else
      // 		    cout << "!! " << " (.<b):  \t  Ju="  << Ju
      // 			 << " !> " << x << "=b" << endl;
      // 		  break;
      // 		}
      // 	      case Bound::BOUND_DOUBLE:
      // 		{
      // 		  double xi = b[i][r].getBound( Bound::BOUND_INF );
      // 		  double xs = b[i][r].getBound( Bound::BOUND_SUP );
      // 		  if( (xi<Ju+EPSILON)&&(Ju<xs+EPSILON) )
      // 		    cout << "    (b<.<b):  \t  binf="<<xi<<" < Ju="  << Ju
      // 			 << " < " << xs << "=bsup" << endl;
      // 		  else
      // 		    cout << "!!  (b<.<b):  \t  binf="<<xi<<" !< Ju="  << Ju
      // 			 << " !< " << xs << "=bsup" << endl;
      // 		  break;
      // 		}
      // 	      }


      // 	  }
      // }


      // sotDEBUG(1) << "Check bounds of " << i << "."<<endl;
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
      } */



  if(! DummyActiveSet::explore(J,b,solution) )
    {      exit(-1);    }

  // DummyActiveSet::detailActiveSet( J,b,solution,2720,0 );
  // DummyActiveSet::detailActiveSet( J,b,solution,2721,1 );
  // DummyActiveSet::detailActiveSet( J,b,solution,2736,0 );
}
