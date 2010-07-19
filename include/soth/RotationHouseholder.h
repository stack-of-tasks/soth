
namespace soth
{

  /* --- ROT HH ------------------------------------------------------------- */
  /* --- ROT HH ------------------------------------------------------------- */
  /* --- ROT HH ------------------------------------------------------------- */

  class RotationHouseHolder
  {
  protected:
    VectorXd h;
    double factor;

  public:
    RotationHouseHolder( void )

    // Init from column j of QR-issued matrix m.
    template< typename MatrixType >
    void init( MatrixType& m, unsigned int j );

    /* --- operator -- */
    template <typename MatrixType >
    MatrixType & multiplyLeft( MatrixType & m );
    template <typename MatrixType >
    MatrixType & multiplyRight( MatrixType & m );

  };

  typedef std::list< RotationGiven > RotationGiven_list_t;
  typedef std::list< RotationHouseHolder > RotationHouseHolder_list_t;

  template< typename MatrixGen >
  void initFromQR( RotationHouseHolder_list_t& hh, MatrixGen& A )
  {
    const unsigned int colA=A.size2();
    const unsigned int rowA=A.size1();
    hh.resize(colA);
    RotationHouseHolder_list_t::iterator iterHH=hh.begin();
    for( unsigned int j=0;j<colA;++j,++iterHH )
      {
	iterHH->init(A,j);
	for( unsigned int i=j+1;i<rowA; ++i ) { A(i,j) = 0; }
      }
  }

  template< typename MatrixGen >
  void prod( const RotationHouseHolder_list_t& hh, MatrixGen& A );
  template< typename MatrixGen >
  void prod( MatrixGen& A,const RotationHouseHolder_list_t& hh );




}; // namespace soth
