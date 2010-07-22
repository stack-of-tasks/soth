#include "soth/Givens.h"

namespace soth
{
  Givens::Givens()
    :i(0), j(0)
  {
  }

  Givens::Givens(double a, double b, int i, int j, double* z)
    :i(i), j(j)
  {
    G.makeGivens(a,b,z);
  }

  void Givens::makeGivens(double a, double b, int i, int j, double* z)
  {
    G.makeGivens(a,b,z);
    this->i = i;
    this->j = j;
  }
}