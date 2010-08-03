#ifndef __SOTH_HCOD__
#define __SOTH_HCOD__

#include "soth/Stage.hpp"
#include "soth/BaseY.hpp"
#include <boost/smart_ptr.hpp>

namespace soth
{


  class HCOD
  {
  protected:
    typedef boost::shared_ptr<soth::Stage> stage_ptr_t;
    typedef std::vector<stage_ptr_t> stage_sequence_t;
    typedef stage_sequence_t::iterator stage_iter_t;
    typedef stage_sequence_t::reverse_iterator stage_riter_t;
    typedef std::vector<VectorXi> activeset_sequence_t;

  public:
    HCOD( unsigned int sizeProblem, unsigned int nbStage = 0 );

    void pushBackStage( const MatrixXd & J, const bound_vector_t & bounds );
    void pushBackStage( const MatrixXd & J, const bound_vector_t & bounds,const VectorXi& Ir0 );
    void pushBackStages( const std::vector<MatrixXd> & J,
			 const std::vector<bound_vector_t> & bounds );

    Stage& stage( unsigned int i );
    const Stage& stage( unsigned int i ) const;
    Stage& operator[] ( unsigned int i ) { return stage(i); }
    const Stage& operator[] ( unsigned int i ) const { return stage(i); }

    void setInitialActiveSet( const VectorXi& Ir0,unsigned int i );
    const VectorXi& getInitialActiveSet( unsigned int i );

    //sizes
    int sizeA() const;
    int rank() const;
    int nbStages() const { return stages.size(); }

    /* --- Decomposition --- */
  public:
    void reset( void );
    void initialize( void );
    void update( const unsigned int & stageUp,const Stage::ConstraintRef & cst );
    void update( stage_iter_t stageIter,const Stage::ConstraintRef & cst );
    void downdate( const unsigned int & stageDown, const unsigned int & row );
    void downdate( stage_iter_t stageIter,const unsigned int & row );
  protected:
    void updateY( const GivensSequence& Yup );

    /* --- Computations --- */
  public:
    void computeSolution( bool compute_u = true );
    void computeLagrangeMultipliers( void );
    double computeStepAndUpdate( void );
    double computeStep( void );
    bool searchAndDowndate( void );
    bool search( void );

    void makeStep( double tau, bool compute_u = true );
    void makeStep( bool compute_u = true );



    /* --- The big one --- */
  public:
    //template< typename VectorGen >
    void activeSearch( VectorXd & u );

    /* --- Tests --- */
  public:
    void show( std::ostream& os, bool check=false );
    bool testRecomposition( std::ostream* os );
    bool testLagrangeMultipliers( std::ostream* os ) const;
    bool testLagrangeMultipliers( std::ostream& os ) const
    { testLagrangeMultipliers(&os); }

    void setNameByOrder( const std::string root = ""  );

  protected:
    HCOD( void ) : Y(0) {};

  protected:public://DEBUG
    unsigned int sizeProblem;
    soth::BaseY Y;
    stage_sequence_t stages;
    activeset_sequence_t initialActiveSets;
    VectorXd solution;

    VectorXd du,Ytu,Ytdu,rho;
    bool isReset,isInit,isSolutionCpt;
  };



}; // namespace soth


#endif // #ifndef __SOTH_HCOD__
