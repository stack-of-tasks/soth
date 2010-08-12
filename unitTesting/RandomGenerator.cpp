#include "RandomGenerator.h"
#include "MatrixRnd.h"
#include "soth/debug.h"
#include <fstream>

namespace soth
{

  using std::endl;

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
	    //if( s==0 ) btype= randu(2,4); // DEBUG
	    //if( s!=0 ) btype= 1; // DEBUG
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
    nc = Random::rand<int>() % 30 +4 ; //%  50 + 6;
    nbStage = randu(1,1+nc/4);

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

	// DEBUG
	if(rankfree[i]==0 ) rankfree[i]++;

      }

  }



  void randomProblem( std::vector<Eigen::MatrixXd> &J,
		      std::vector<soth::bound_vector_t> &b,
		      bool verbose )
  {
    int NB_STAGE,NC;
    std::vector<int> NR,RANKLINKED,RANKFREE;
    generateRandomProfile(NB_STAGE,RANKFREE,RANKLINKED,NR,NC);

    /* Initialize J and b. */
    generateDeficientDataSet(J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC);



    if( verbose )
      {
	for( unsigned int i=0;i<NB_STAGE;++i )
	  { std::cout << RANKFREE[i]<<"/"<<NR[i]<<",  "; }

	std::cout << "NB_STAGE=" << NB_STAGE <<",  NC=" << NC << endl;
	std::cout << endl;
	for( unsigned int i=0;i<NB_STAGE;++i )
	  {
	    std::cout << "J"<<i+1<<" = " << (soth::MATLAB)J[i] << std::endl;
	    std::cout << "e"<<i+1<< " = " << b[i] << ";"<<std::endl;
	  }
      }
  }




  void readProblemFromFile( const std::string name,
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

  void readProblemFromFile( const std::string name,
			    std::vector<Eigen::MatrixXd> &J,
			    std::vector<soth::bound_vector_t> &b )
  {
    int NB_STAGE;
    std::vector<int> NR;
    int NC;
    readProblemFromFile( name,J,b,NB_STAGE,NR,NC);
  }

  void writeProblemToFile( const std::string name,
			   const std::vector<Eigen::MatrixXd> &J,
			   const std::vector<soth::bound_vector_t> &b,
			   const int& NB_STAGE,
			   const std::vector<int> & NR,
			   const int& NC )

  {
    std::ofstream fout(name.c_str());
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

  void writeProblemToFile( const std::string name,
			   const std::vector<Eigen::MatrixXd> &J,
			   const std::vector<soth::bound_vector_t> &b )
  {
    int NB_STAGE;
    std::vector<int> NR;
    int NC;
    writeProblemToFile( name,J,b,NB_STAGE,NR,NC);
  }

}
