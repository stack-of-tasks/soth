namespace soth
{

  /* --- ROT GR ------------------------------------------------------------- */
  /* --- ROT GR ------------------------------------------------------------- */
  /* --- ROT GR ------------------------------------------------------------- */

  class RotationGiven
  {
  protected:
    double c,s;
    unsigned int i,j;

  public:
    RotationGiven( double a, double b, unsigned int i, unsigned int j );
    RotationGiven( void );

    template< typename MatrixType >
    void init( const MatrixType& m, unsigned int i1, unsigned int i2, unsigned int j );
    template< typename MatrixType >
    void init( unsigned int i, unsigned int j1, unsigned int j2, const MatrixType& m );


    /* --- operator -- */
    template <typename MatrixType >
    MatrixType & multiplyLeft( MatrixType & m );
    template <typename MatrixType >
    MatrixType & multiplyRight( MatrixType & m );

  };



}; // namespace soth
