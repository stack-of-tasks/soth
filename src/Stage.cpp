#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>

namespace bnu = boost::numeric::ublas;
typedef bnu::Matrix<double> MatrixXd;
typedef bnu::Vector<double> VectorXd;
typedef bnu::index_matrix<> Indirect;

typedef bnu::matrix_indirect< MatrixXd,Indirect > SubMatrixXd;
typedef bnu::triangular_adaptator< SubMatrixXd > TriSubMatrixXd;

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


  /* --- STAGE -------------------------------------------------------------- */
  /* --- STAGE -------------------------------------------------------------- */
  /* --- STAGE -------------------------------------------------------------- */

  class Stage
  {
  protected:
    unsigned int nr,nc; // nr=nbCols, nc=nbRows.

    const MatrixXd & J_;
    MatrixXd W_;
    MatrixXd ML_;

    SubMatrixXd M,L;
    SubMatrixXd L0sq;
    TriSubMatrixXd L0; // L0 = tri(L0sq)

    TriMatrixXd Ldamp;

    // fullRankRows = Ir. defRankRows = In.
    Indirect& Ir,In; // Ir = L0sq.indirect1() -- In = 
    Indirect unactiveRows;

    unsigned int sizeM,sizeL; // sizeL = card(Ir)

    // W = W_( :,[In Ir] ).
    // M = ML_( [In Ir],0:sizeM-1 ).
    // L = ML_( Ir,sizeM:sizeM+sizeL-1  ).
    // Lo = [0;L] = ML_( [In Ir],sizeM:sizeM+sizeL-1  ).

    // W_*ML_ = W*[M [ zeros(In.size(),rank); L ] ]


  public:
    HStage( const MatrixXd & J );

    /* --- INIT ------------------------------------------------------------- */

    void initialCOD( const unsigned int previousRank,
		     const Indirect & initialIr,
		     MatrixXd & Y);
    /*
      ML=J(initIr,:)*Y;
      rank=Ir.size();  Ir=1:rank;
      M=submatrix(ML,1:previousRank); L=submatrix(ML,previousRank+1:end);

      A=columnMajor(L)  // A==L'
      qr(A);
      RotationHouseHolder_list_t Yup( A );
      Y=Y*Yup;

      for i=rank:-1:1
        if( L(i,i)!= 0 ) break;
	nullifyLineDeficient( i );

     */


    void nullifyLineDeficient( const unsigned int r );
    /*
      Jd = L.row(r);
      foreach i in rank:-1:1
        if( Jd(i)==0 ) continue;
	gr= GR(L(Ir(i),i),Jd(i),i,r );
	L=gr*L;
	W=W*gr';
      Ir >> r;
      In << r;
     */

    /* --- DOWN ------------------------------------------------------------- */

    // Return true if the rank re-increase operated at the current stage.
    bool downdate( const unsigned int position,
		   RotationGiven_list_t & Ydown );
    /*
      gr = Neutralize row <position> in W <position>
      L=gr'*L
      bool res;
      if( L(In.last(),0) == 0
        // Rank deficience: regularize Hessenberg
	Ydown = regularizeHessenberg
	res=false;
      else
        // No rank dec: quit
	Ir << In.pop();
	res=true;
      Ir >> r; In >> r;
      Unused << r;
      return res;
     */


    // Return true if the rank decrease operated at the current stage.
    bool propagateDowndate( RotationGiven_list_t & Ydown,
			    bool decreasePreviousRank );
    /*
      M=M*Ydown;
      if(! decreasePreviousRank ) return true;
      L.indices2().push_front( M.indice2().pop_back() );

      foreach i in In
        if( L(i,0) == 0 continue;
	Ir << i; In >> i;
        return true;

      Ydown += regularizeHessenberg
      return false
     */

    void regularizeHessenberg( RotationGiven_list_t & Ydown,
			       unsigned int i0 = 0 );
    /*
      for i=i0:rank-1
        gr = GR( L(Ir(i),i),L(Ir(i),i+1),i,i+1 );
	L = L*GR;
	Ydown.push_back( gr );
     */

    /* --- UPD -------------------------------------------------------------- */

    /* TODO
      if stage(sup).update( Jup,Yup )
        while stage(sup++).propagateUpdate( Yup,false
      else ..

     */

    // Return true if the rank re-decrease operated at the current stage.
    bool update( const VectorXd & Jup,
		 RotationGiven_list_t & Yup );
    /*
      Inew = Unused.pop();
      Row JupY = row(Inew);
      JupU = Jup*Y;
      double norm2=0; double rankJ=0;
      for i=n:-1:1
        norm2+=JupY(i)^2;
	if norm2!=0 rankJ=i; break;

      Ir << Inew
      W(Inew,Inew)=1;
      if rankJ>sizeM+rank
        // Rank increase
        for i=rankJ:-1:sizeM+rank+1
	  Gr gr; gr.init( JupY,i,i-1,0 ); prod(JupY,gr);
	  Yup.push_back( gr );
	return false;
      else
        // No rank increase;
	nullifyLineDeficient(Inew);
	return true;
     */


    // Return true if the rank decrease operated at the current stage.
    bool propagateUpdate( RotationGiven_list_t & Yup,
			  bool decreasePreviousRank );
    /*
      ML=ML*Yup;
      increaseSizeM();
      if( decreasePreviousRank ) return true
      for i=1:rank
        if L(i,i) == 0
	  // Rank re-decrease
	  nullifyLineDeficient( i );
	  regularizeHessenberg( Yup,i+1 );
	  return true;
      return false;
     */


    /* --- SOLVE ------------------------------------------------------------ */

    void solve( const VectorXd & e, VectorXd res );
    void solveTranspose( const VectorXd & e, VectorXd res );

    void damp( const double & damping );

  public:
    /* --- ACCESSORS --- */

    MatrixSubf M();
    MatrixSubf L();
    MatrixSubf Lo();

    const_MatrixSubf M() const ;
    const_MatrixSubf L() const ;
    const_MatrixSubf Lo() const ;

    RowSubf row( const unsigned int r );
    RowSubf rown( const unsigned int r ); // row rank def
    RowSubf rowf( const unsigned int r ); // row rank full

    /* --- MODIFIORS --- */
    void increaseSizeM();
    void decreaseSizeM();

  };




}; // namespace soth
