#ifndef __SOTH_HCOD__
#define __SOTH_HCOD__

#include "soth/Stage.hpp"
#include "soth/BaseY.hpp"
#include <boost/smart_ptr.hpp>

namespace soth
{


  class HCOD
  {
  public:
    HCOD( unsigned int sizeProblem, unsigned int nbStage = 0 );

    void pushBackStage( const MatrixXd & J, const bound_vector_t & bounds );
    void pushBackStage( const MatrixXd & J, const bound_vector_t & bounds,const VectorXi& Ir0 );

    Stage& stage( unsigned int i );
    const Stage& stage( unsigned int i ) const;

    void setInitialActiveSet( const VectorXi& Ir0,unsigned int i );
    const VectorXi& getInitialActiveSet( unsigned int i );

    void updateY( const GivensSequence& Yup );

    //sizes
    int sizeA() const;
    int rank() const;

    void reset( void );
    void solve( void );
    void computeLambda();

    void update( const unsigned int & stageUp,const Stage::ConstraintRef & cst );
    void downdate( const unsigned int & stageDown, const unsigned int & row );

    void show( std::ostream& os, bool check=false );

  template< typename VectorGen >
  void activeSearch( VectorGen & u );


  protected:
    HCOD( void ) : Y(0) {};

  protected:
    typedef boost::shared_ptr<soth::Stage> stage_ptr_t;
    typedef std::vector<stage_ptr_t> stage_sequence_t;
    typedef std::vector<VectorXi> activeset_sequence_t;

  protected:
    unsigned int sizeProblem;
    soth::BaseY Y;
    stage_sequence_t stages;
    activeset_sequence_t initialActiveSets;
    VectorXd solution;
    VectorXd lambda; //lagrange multiplier

  };



}; // namespace soth


#endif // #ifndef __SOTH_HCOD__
