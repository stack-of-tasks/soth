namespace soth
{
  class Random
  {
    static const unsigned int SOTH_RND_MAX = 4278255361U;
    static const unsigned int MULT = 1884103651;
    static unsigned int current;

    static unsigned int next();

  public:
    static void setSeed(unsigned int newSeed);

    template <typename ReturnType>
      static ReturnType rand();

    template <typename ReturnType>
      static ReturnType randMax();
  };

  unsigned int Random::current = 33331;

  unsigned int Random::next()
  {
    static const unsigned long long int m = static_cast<unsigned long long int>(MULT);
    const unsigned long long int c = m * current;
    current = c % SOTH_RND_MAX;
    return current;
  }

  void Random::setSeed(unsigned int newSeed) {current = newSeed;}

  template <> unsigned int Random::rand() {return next();}
  template <> unsigned int Random::randMax() {return SOTH_RND_MAX;}
  template <> int Random::rand() {return next()>>1;}
  template <> int Random::randMax() {return SOTH_RND_MAX>>1;}
  template <> double Random::rand() {return static_cast<double>(next())/SOTH_RND_MAX;}
  template <> double Random::randMax() {return 1.;}


  template<typename Scalar> 
    struct rnd_traits
    {
      static Scalar max() {return Random::randMax<Scalar>();}
      static Scalar default_min() {return Scalar(-1);}
      static Scalar default_max() {return Scalar(1);}
    };

  template<>
    struct rnd_traits<int>
    {
      static int max() {return Random::randMax<int>();}
      static int default_min() {return -10;}
      static int default_max() {return 11;}
    };

  class MatrixRnd
  {
  public:
    template<typename Derived>
      static MatrixBase<Derived>&
      randomize(MatrixBase<Derived>& m,
		typename Derived::Scalar min= rnd_traits<typename Derived::Scalar>::default_min(),
		typename Derived::Scalar max= rnd_traits<typename Derived::Scalar>::default_max())
      {
	for ( typename Derived::Index i=0; i<m.rows(); ++i)
	  {
	    for ( typename Derived::Index j=0; j<m.cols(); ++j)
	      m(i,j) = static_cast<typename Derived::Scalar>((max-min)* Random::rand<double>()) + min;
	  }
	return m;
      }
  };

}; // namespace soth.
