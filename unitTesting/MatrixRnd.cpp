#include "MatrixRnd.hpp"

#ifdef WIN32
inline double round( double d )
{	return floor( d + 0.5 );	}
#endif /* WIN32 */

namespace soth
{
  unsigned int Random::current = 33331;

  unsigned int Random::next()
  {
    static const unsigned long long int m = static_cast<unsigned long long int>(MULT);
    const unsigned long long int c = m * current;
    current = c % SOTH_RND_MAX;
    return current;
  }

  void Random::setSeed(unsigned int newSeed) {current = newSeed;}


  // Simulate a white noise with mean 0 and var 1.
  double whiteNoise(void)
  {
    const int ACC = 100;
    double x=0;
    for( int i=0;i<ACC;++i ) x=x+Random::rand<double>();
    return (x-ACC/2.)*sqrt(12.0/ACC);
  }
  // Simulate any white noise.
  int whiteNoise( int mean,double var )
  {
    double x=whiteNoise()*var+mean;
    return std::max(0,(int)round(x));
  }
  // Simulate a discrete uniform law inside [ bmin,bmax ] (each bound
  // being reached).
  int randu( int bmin,int bmax )
  {
    assert( bmin<bmax );
    return floor((bmax-bmin+1)*Random::rand<double>()+bmin);
  }


}
