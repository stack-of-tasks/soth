/*
 *  Copyright
 */

#include "soth/Stage.hpp"
#include <iostream>
#include <Eigen/QR>
#include <Eigen/Jacobi>
#include <boost/smart_ptr.hpp>


int main (int argc, char** argv)
{
  // const int NB_STAGE = 5;
  // const int RANK[] = { 3, 5, 3, 5, 3 };
  // const int NR[] = { 5, 5, 5, 5, 8 };
  // const int NC = 20;

  const unsigned int NB_STAGE = 3;
  const unsigned int RANK[] = { 3, 4, 3, 5, 3 };
  const unsigned int NR[] = { 5, 4, 5, 5, 8 };
  const unsigned int NC = 12;
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
      for( unsigned int j=0;j<NR[i];++j ) b[i][j] = j*i*0.5;

      std::cout << "J"<<i<<" = " << (soth::MATLAB)J[i] << std::endl;
    }

  /* SOTH structure construction. */
  soth::BaseY Y(NC);
  typedef boost::shared_ptr<soth::Stage> stage_ptr_t;
  typedef std::vector<stage_ptr_t> stage_list_t;
  stage_list_t stages(NB_STAGE);
  for( unsigned int i=0;i<NB_STAGE;++i )
    {
      std::cout <<" --- STAGE " <<i<< " --------------------------------- " << std::endl;
      /* Compute the initial COD for each stage. */
      stages[i] = stage_ptr_t(new soth::Stage( J[i],b[i],Y ));
      stages[i]->computeInitialCOD(soth::Stage::allRows(),Y);
      Eigen::MatrixXd Jrec; stages[i]->recompose(Jrec);
      std::cout << "Jrec" <<i<<" = " << (soth::MATLAB)Jrec << std::endl;
    }

  Eigen::MatrixXd rec;
  //stages[0]->recompose(rec);
}
