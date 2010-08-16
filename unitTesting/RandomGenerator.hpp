
#ifndef __SOTH_RANDOM_GENERATOR__
#define __SOTH_RANDOM_GENERATOR__

#include "soth/Algebra.hpp"
#include "soth/Bound.hpp"
#include <vector>

namespace soth
{
  /* Being specified the size of the problem NC, the number of stage NB_STAGE,
   * the size of each size NR, the rank of the stage wrt. the previous stages
   * RANKFREE = rank(Ji.Pi-1), and the rank due to the previous stage
   * RANKLINKED = rank(Ji) - RANKFREE, the function generates in J and b
   * a random problems with good conditionning numbers (ie =0 XOR >>0).
   */
  void generateDeficientDataSet( std::vector<Eigen::MatrixXd> &J,
				 std::vector<soth::bound_vector_t> &b,
				 const int NB_STAGE,
				 const std::vector<int> & RANKFREE,
				 const std::vector<int> & RANKLINKED,
				 const std::vector<int> & NR,
				 const int NC );

  /* Generated randomly a profile of problem, ie sizes and ranks. The ouput
   * of this functions are to be sent to the previous function. */
  void generateRandomProfile(int & nbStage,
			     std::vector<int>& rankfree,
			     std::vector<int>& ranklinked,
			     std::vector<int>& nr,
			     int & nc );

  void randomProblem( std::vector<Eigen::MatrixXd> &J,
		      std::vector<soth::bound_vector_t> &b );


  /* --- IN/OUT PROBLEMS --- */
  void readProblemFromFile( const std::string name,
			    std::vector<Eigen::MatrixXd> &J,
			    std::vector<soth::bound_vector_t> &b,
			    int& NB_STAGE,
			    std::vector<int> & NR,
			    int& NC );
  void readProblemFromFile( const std::string name,
			    std::vector<Eigen::MatrixXd> &J,
			    std::vector<soth::bound_vector_t> &b );

  void writeProblemToFile( const std::string name,
			   const std::vector<Eigen::MatrixXd> &J,
			   const std::vector<soth::bound_vector_t> &b,
			   const int& NB_STAGE,
			   const std::vector<int> & NR,
			   const int& NC );

  void writeProblemToFile( const std::string name,
			   const std::vector<Eigen::MatrixXd> &J,
			   const std::vector<soth::bound_vector_t> &b );
};

#endif // #ifndef __SOTH_RANDOM_GENERATOR__
