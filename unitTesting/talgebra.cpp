/*
 *  Copyright
 */

#include "soth/Algebra.h"
#include <iostream>

bool soth::MATLAB::fullPrec = true;

int main (int argc, char** argv)
{
  soth::MatrixXd A(4,5);
  soth::randMatrix(A,9,5);
  std::cout << (soth::MATLAB)A << std::endl;



}
