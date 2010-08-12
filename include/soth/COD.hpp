#ifndef __SOTH_COD__
#define __SOTH_COD__

#include "soth/Algebra.h"
#include <Eigen/SVD>

namespace soth
{

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

    void disp( bool check=false, double EPSILON=1e-4 )
    {
      if( rank==0 )
	{
	  sotDEBUG(5) << "ULV: Empty rank."  << std::endl;
	  return;
	}

      sotDEBUG(5) << "U = " << (MATLAB)(MatrixXd)qru.householderQ() << std::endl;
      sotDEBUG(5) << "R = " << (MATLAB)(MatrixXd)Rt.transpose() << std::endl;
      sotDEBUG(5) << "V = " << (MATLAB)(MatrixXd)qrv.householderQ() << std::endl;
      sotDEBUG(5) << "L = " << (MATLAB)(MatrixXd)qrv.matrixQR().topRows(rank).triangularView<Upper>().transpose() << std::endl;
      sotDEBUG(5) << "Pcv = " << (MATLAB)(MatrixXd)(qru.colsPermutation()) << std::endl;
      sotDEBUG(5) << "Pcu = " << (MATLAB)(MatrixXd)(qrv.colsPermutation()) << std::endl;

      if( check )
	{
	  MatrixXd L0(NR,NC); L0.setZero();
	  L0.topLeftCorner(rank,rank)
	    = qrv.matrixQR().topRows(rank).triangularView<Upper>().transpose();
	  MatrixXd Arec = matrixU() * L0 * matrixV().transpose();
	  sotDEBUG(5) << "Arec = " << (MATLAB)Arec << std::endl;
	  assert( (Ainit-Arec).norm()< EPSILON );
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
      sotDEBUG(5) << "VLUe = " << (MATLAB)solv << std::endl;
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


}


#endif // #ifndef __SOTH_COD__
