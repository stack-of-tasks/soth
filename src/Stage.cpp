#include "soth/Stage.hpp"

namespace soth
{

  Stage::
  Stage( const MatrixXd & inJ, const bound_vector_t & inbounds,BaseY & inY )
    : J(inJ), bounds(inbounds)
    ,Y(inY)
    ,nr(J.rows()),nc(J.cols())
    ,W_(nr,nr),ML_(nr,nc),e_(nr)
    ,sizeM(0),sizeL(0)
  {
    assert( bounds.size() == J.rows() );
  }


  void Stage::
  computeInitialCOD( const unsigned int previousRank,
		     const Indirect & initialIr )
  {



  }




}; // namespace soth
