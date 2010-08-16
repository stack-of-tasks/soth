#include <Eigen/Core>
#include "soth/SubMatrix.hpp"
#include <iostream>
#include <string>

using namespace Eigen;

void testPermutedBasic()
{
  MatrixXi m(4,5);
  Map<MatrixXi> map(m.data(), m.size(), 1);
  map = VectorXi::LinSpaced(m.size(), 0, m.size()-1);
  std::cout << "original matrix:" << std::endl << std::endl;
  std::cout << m << std::endl;

  SubMatrix<MatrixXi> p(m,true,true);
  std::cout << "creating a permuted matrix with default permutation" << std::endl;
  std::cout << p << std::endl << std::endl;

  std::cout << "permuting row 1 and 2" << std::endl;
  p.permuteRows(1,2);
  std::cout << p << std::endl << std::endl;
  assert( p(1,1)==6 );


  std::cout << "permuting col 2 and 4 then 1 and 4" << std::endl;
  p.permuteCols(2,4);
  p.permuteCols(1,4);
  std::cout << "col permutation indices are " << p.getColIndices().transpose() << " and the permuted matrix is" << std::endl;
  std::cout << p << std::endl << std::endl;
  assert( p(1,1)==10 );  assert( p(1,4)==6 );

  std::cout << "*****Read access******" << std::endl;
  std::cout << "coef (1,4) is " << p(1,4) << std::endl;
  std::cout << "top right corner of size (2,2) is " << std::endl;
  std::cout << p.topRightCorner(2,2) << std::endl;
  assert(p.topRightCorner(2,2)(1,0)==14); assert(p.topRightCorner(2,2)(0,1)==4);

  std::cout << "and a triangular lower view" << std::endl;
  MatrixXi pt = p.triangularView<Lower>();
  std::cout << pt << std::endl;
  assert(pt(3,3)==15);  assert(pt(3,4)==0);

  std::cout << "*****Write access******" << std::endl;
  std::cout << "set coef (2,1) to -1 and col(2) to -2" << std::endl;
  p(2,1) = -1;
  p.col(2).setConstant(-2);
  assert( m(1,2)==-1 );
  assert( m(0,4)==-2 );  assert( m(3,4)==-2 );

  std::cout << "permuted matrix:" << std::endl;
  std::cout << p << std::endl;
  std::cout << "original matrix:" << std::endl;
  std::cout << m << std::endl;

  std::cout << std::endl << std::endl;
  std::cout << "----Now a row permuted only----" << std::endl;
  SubMatrix<MatrixXi, RowPermutation> pr(m,true);
  std::cout << "creating a row permuted matrix with default permutation" << std::endl;
  std::cout << pr << std::endl << std::endl;

  std::cout << "permuting row 1 and 2" << std::endl;
  assert( pr(1,1)==5 );
  pr.permuteRows(1,2);
  std::cout << pr << std::endl << std::endl;
  assert( pr(1,1)==6 );

  std::cout << "*****Read access******" << std::endl;
  std::cout << "coef (1,4) is " << pr(1,4) << std::endl;
  assert( pr(1,4)==-2 );
  std::cout << "top left corner of size (2,2) is " << std::endl;
  std::cout << pr.topLeftCorner(2,2) << std::endl;
  assert(pr.topLeftCorner(2,2)(1,0)==2); assert(pr.topLeftCorner(2,2)(0,1)==4);

  std::cout << "and a triangular lower view" << std::endl;
  pt = pr.triangularView<Lower>();
  std::cout << pt << std::endl;
  assert(pt(3,3)==15);  assert(pt(3,4)==0);

  std::cout << std::endl << std::endl;
  std::cout << "----Now a col submatrix----" << std::endl;
  SubMatrix<MatrixXi, ColPermutation>::ColIndices ind(3);
  ind << 2,4,1;
  SubMatrix<MatrixXi, ColPermutation> pc(m);
  pc.setColIndices(ind);
  std::cout << "original matrix:" << std::endl << std::endl;
  std::cout << m << std::endl;
  std::cout << "sub matrix" << std::endl;
  std::cout << pc << std::endl << std::endl;
  assert(pc.cols()==3);         assert( pc(0,0) == m(0,2) );
  assert( pc(0,1) == m(0,4) );  assert( pc(0,2) == m(0,1) );

}

void testSubVector()
{

  VectorXi v = VectorXi::LinSpaced(8,0,7);
  std::cout << "original vector:" << std::endl;
  std::cout << v.transpose() << std::endl << std::endl;

  VectorXi row(4); row << 2,3,5,0;
  SubMatrix<VectorXi, RowPermutation> pv(v);
  pv.setRowIndices(row);
  std::cout << "subVector:" << std::endl;
  std::cout << pv.transpose() << std::endl << std::endl;
  assert(pv[0]==2);  assert(pv[3]==0);

  std::cout << "*****Read access******" << std::endl;
  std::cout << "coef 2 is " << pv[2] << std::endl;
  assert(pv[2]==5);

  std::cout << "*****Operation******" << std::endl;
  std::cout << "subVector:" << std::endl;
  VectorXi res = pv.transpose() + Vector4i::Ones().transpose();
  std::cout << res << std::endl << std::endl;
  assert( res[1] == 4 );

  VectorXi row2(4); row2 << 7,6,1,4;
  SubMatrix<VectorXi, RowPermutation> pv2(v);
  pv2.setRowIndices(row2);
  pv2 = pv;
  std::cout << pv2.transpose() << std::endl;
  std::cout << v.transpose() << std::endl;

  pv.popRowFront();
  std::cout << pv.transpose() << std::endl;
  pv.pushRowFront(6);
  std::cout << pv.transpose() << std::endl;
  assert( pv[0] == 3 );  assert( pv[1] == 3 );
}


void testDoublePerm()
{
  std::cout << "\n\n **** \n **** Proxy\n **** "<< std::endl;

  MatrixXi m(4,5);
  Map<MatrixXi> map(m.data(), m.size(), 1);
  map = VectorXi::LinSpaced(m.size(), 0, m.size()-1);
  std::cout << "original matrix:" << std::endl << std::endl;
  std::cout << m << std::endl;

  VectorXi idx(3); idx << 3,1,2;
  SubMatrix<MatrixXi> p(m,&idx,&idx);
  std::cout << "double permuted matrix:" << std::endl << "."<<p<<"."<<std::endl;
  assert( p(1,2) == 9 );

  p.pushRowFront( p.popRowBack() );
  std::cout << "pass last in first:" << std::endl << p<<std::endl;
  std::cout << "ir = [" << p.getRowIndices().transpose() << " ] " << std::endl;
  std::cout << "ic = [" << p.getColIndices().transpose() << " ] " <<  std::endl;
  assert( &p.getRowIndices()==&p.getColIndices() );
  assert( p(1,1) == 15 ); assert( p(0,2) == 6 );
}


void speedTest()
{
  std::cout << "\n\n **** \n **** Speed \n **** "<< std::endl;
  int N=1000;
  int n[] = {5,10,50,100,250};
  for (int i=0; i<5; ++i)
  {
    std::cout << "size " << n[i] << std::endl;
    MatrixXd A = MatrixXd::Random(n[i],n[i]);
    MatrixXd B = MatrixXd::Random(n[i],n[i]);
    MatrixXd C1(n[i],n[i]);
    MatrixXd C2(n[i],n[i]);
    MatrixXd C3(n[i],n[i]);
    MatrixXd C4(n[i],n[i]);
    SubMatrix<MatrixXd> Prc(A,true,true);
    SubMatrix<MatrixXd, RowPermutation> Pr(A,true);
    SubMatrix<MatrixXd, ColPermutation> Pc(A,true);

    double dummy;
    clock_t start, stop;
    double total;

    dummy = 0;
    start = clock();
    for (int j=0; j<N; ++j)
    {
      C1.noalias() = A*B;
      dummy += C1(0,0);
    }
    stop = clock();
    total = static_cast<double>(stop-start)/(CLOCKS_PER_SEC*N)*1000;
    std::cout << "normal mult    : " << total  << "ms" << std::endl;

    dummy = 0;
    start = clock();
    for (int j=0; j<N; ++j)
    {
      C2.noalias() = Prc*B;
      dummy += C2(0,0);
    }
    stop = clock();
    total = static_cast<double>(stop-start)/(CLOCKS_PER_SEC*N)*1000;
    std::cout << "Perm*Normal    : " << total  << "ms" << std::endl;

    dummy = 0;
    start = clock();
    for (int j=0; j<N; ++j)
    {
      C3.noalias() = Pr*B;
      dummy += C3(0,0);
    }
    stop = clock();
    total = static_cast<double>(stop-start)/(CLOCKS_PER_SEC*N)*1000;
    std::cout << "RowPerm*Normal : " << total  << "ms" << std::endl;

    dummy = 0;
    start = clock();
    for (int j=0; j<N; ++j)
    {
      C4.noalias() = Pc*B;
      dummy += C4(0,0);
    }
    stop = clock();
    total = static_cast<double>(stop-start)/(CLOCKS_PER_SEC*N)*1000;
    std::cout << "ColPerm*Normal : " << total  << "ms" << std::endl;

    std::cout << std::endl;
  }
}

void testCol()
{
  VectorXi row(3); row << 2,3,0;
  MatrixXi m(4,5);
  Map<MatrixXi> map(m.data(), m.size(), 1);
  map = VectorXi::LinSpaced(m.size(), 0, m.size()-1);

  MatrixXi::ColXpr m1 = m.col(1);
  SubMatrix< MatrixXi::ColXpr,RowPermutation > m1i( m1,row );


  SubCol<MatrixXi> l(m,row,1);
  std::cout << l << std::endl;
  assert( (l-m1i).sum() == 0 );

}


int main(int argc, char**argv)
{
  //testPermutedBasic();
  //testSubVector();
  //testDoublePerm();
  testCol();
  if( (argc>1)&&(std::string(argv[1]) == "-speed") ) speedTest();

  std::cout << "\n\n ... Everything is working fine.\n\n\n" << std::endl;
}
