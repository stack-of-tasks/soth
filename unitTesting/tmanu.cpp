/*
 *  Copyright
 */
#define SOTH_DEBUG
#define SOTH_DEBUG_MODE 45
#include "soth/debug.hpp"
#include "soth/HCOD.hpp"
#include "soth/COD.hpp"
#include <sys/time.h>
#include <iostream>
#include <vector>

using namespace soth;
using std::endl;
using std::cout;
using std::cerr;

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */


#include <boost/program_options.hpp>

bool checkColumn( const Eigen::MatrixXd& A,
		  const Eigen::MatrixXd& B,
		  const int index )
{
  const int NC = A.rows(),NA=A.cols(),NB=B.cols(), NX = NA+NB;

  std::vector<Eigen::MatrixXd> J(3);
  std::vector<soth::VectorBound> b(3);

  /* SOTH structure construction. */
  soth::HCOD hcod(NX,3);

  J[0].resize(NA,NX); b[0].resize(NA);
  J[0].leftCols(NA).setIdentity();
  J[0].rightCols(NB).fill(0);
  b[0].fill( soth::Bound(0,soth::Bound::BOUND_INF) );
  hcod.pushBackStage( J[0],b[0] );

  J[1].resize(NC,NA+NB); b[1].resize(NC);
  J[1].leftCols(NA) = A;  J[1].rightCols(NB) = B;
  b[1].fill(0);
  hcod.pushBackStage( J[1],b[1] );

  J[2].resize(1,NA+NB); b[2].resize(1);
  J[2].fill(0); J[2](0,index)=1; b[2].fill(0);

  hcod.setNameByOrder("stage_");

  const double dampingFactor = 0.0;
  hcod.setDamping(dampingFactor);
  hcod.stage(0).damping(0);
  hcod.setInitialActiveSet();

  VectorXd solution(NX);
  hcod.activeSearch( solution );
  if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );
  cout << "x = " << (MATLAB)solution << endl;
  cout << "actset = "; hcod.showActiveSet(std::cout);

  cout << "res = " << solution[index] << endl;

  return false;
}



int main (int argc, char** argv)
{
# ifndef NDEBUG
  sotDebugTrace::openFile();
#endif

  const int NB_STAGE =2,NC = 6;
  const int NA=10,NB=10;

  std::vector<Eigen::MatrixXd> J(NB_STAGE);
  std::vector<soth::VectorBound> b(NB_STAGE);;

  Eigen::MatrixXd A = MatrixXd::Random(NC,NA);
  Eigen::MatrixXd B = MatrixXd::Random(NC,NB);
  Eigen::VectorXd Aa = VectorXd::Random(NC);

  /* SOTH structure construction. */
  soth::HCOD hcod(NA+NB,2);

  J[0].resize(NA,NA+NB); b[0].resize(NA);
  J[0].leftCols(NA).setIdentity();
  J[0].rightCols(NB).fill(0);
  b[0].fill( soth::Bound(0,soth::Bound::BOUND_INF) );
  hcod.pushBackStage( J[0],b[0] );

  J[1].resize(NC,NA+NB); b[1].resize(NC);
  J[1].leftCols(NA) = A;  J[1].rightCols(NB) = B;
  for( int i=0;i<NC;++i) b[1][i] = Aa[i];
  hcod.pushBackStage( J[1],b[1] );

  hcod.setNameByOrder("stage_");

  const double dampingFactor = 0.0;
  hcod.setDamping(dampingFactor);
  hcod.stage(0).damping(0);
  hcod.setInitialActiveSet();

  cout << "A = " << (MATLAB)A << endl;
  cout << "Aa = " << (MATLAB)Aa << endl;
  cout << "B = " << (MATLAB)B << endl;

  VectorXd solution(NA+NB);
  hcod.activeSearch( solution );
  if( sotDEBUGFLOW.outputbuffer.good() ) hcod.show( sotDEBUGFLOW.outputbuffer );
  cout << "x = " << (MATLAB)solution << endl;
  cout << "actset = "; hcod.showActiveSet(std::cout);

  cout << "res = " << (J[1]*solution - Aa).norm() << endl;


}
