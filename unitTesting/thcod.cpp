/*
 *  Copyright
 */

#include "soth/HCOD.hpp"
#include "soth/debug.h"


int main (int argc, char** argv)
{
  sotDebugTrace::openFile();
  const int NB_STAGE = 3;
  const int RANK[] = { 1, 4, 3, 5, 3 };
  const int NR[] = { 1, 4, 5, 5, 8 };
  const int NC = 12;

  /* Initialize J and b. */
  std::vector<Eigen::MatrixXd> J(NB_STAGE);
  std::vector<soth::bound_vector_t> b(NB_STAGE);
  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      Eigen::MatrixXd Xhi,Jfr;
      soth::randMatrix(Xhi,NR[i],RANK[i]);
      soth::randMatrix(Jfr,RANK[i],NC);
      J[i]=Xhi*Jfr;
      b[i].resize(NR[i]);
      for( unsigned int j=0;j<NR[i];++j ) b[i][j] = j*(i+1)*0.5;

      /* Introduce a deficience of J1 du to J0. */
      //if( i==0 ) J[0].row(1) = Eigen::MatrixXd::Random(1,NC);
      if( i==0 ) J[0].row(NR[0]-1) = Eigen::MatrixXd::Random(1,NC);
      if( i==1 ) J[1].row(2) = Eigen::MatrixXd::Random(1,NR[0])*J[0];
      if( i==1 ) J[1].row(3) = Eigen::MatrixXd::Random(1,NR[0])*J[0];

      std::cout << "J"<<i<<" = " << (soth::MATLAB)J[i] << std::endl;
      std::cout << "e"<<i;
      for( unsigned int j=0;j<NR[i];++j )
	std::cout << b[i][j].getBound(soth::Bound::BOUND_TWIN) << "   ";
      std::cout << std::endl;
    }
  //b[0][NR[0]-1] = soth::Bound(-0.5, soth::Bound::BOUND_INF);

  /* SOTH structure construction. */
  soth::HCOD hcod(NC,NB_STAGE);
  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      hcod.pushBackStage( J[i],b[i] );
      assert(NR[i]>0);
      if (NR[i]>1)
        hcod.setInitialActiveSet( Eigen::VectorXi::LinSpaced(0, NR[i]-1, NR[i]),i);
      else
        hcod.setInitialActiveSet( Eigen::VectorXi::Zero(1), i);
    }

  hcod.solve();
  //hcod.show(std::cout,true);
}
