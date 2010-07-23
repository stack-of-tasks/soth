/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Copyright Projet Lagadic, 2005
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *
 * File:      sotDebug.h
 * Project:   STack of Tasks
 * Author:    Nicolas Mansard
 *
 * Version control
 * ===============
 *
 *  $Id$
 *
 * Description
 * ============
 *
 * Macro de trace et de debugage
 *
 *   - TRACAGE:    TRACE et ERROR_TRACE fonctionnent comme des printf
 * avec retour chariot en fin de fonction.
 *                 CERROR et CTRACE fonctionnent comme les flux de sortie
 * C++ cout et cerr.
 *   - DEBUGAGE:   DEBUG_TRACE(niv,  et DERROR_TRACE(niv, fonctionnent
 * comme des printf, n'imprimant que si le niveau de trace 'niv' est
 * superieur au mode de debugage SOTH_DEBUG_MODE.
 *                 CDEBUG(niv) fonctionne comme le flux de sortie C++ cout.
 *                 DEBUG_ENABLE(niv) vaut 1 ssi le niveau de tracage 'niv'
 * est superieur au  mode de debugage DEBUG_MODE. Il vaut 0 sinon.
 *   - PROG DEFENSIVE: DEFENSIF(a) vaut a ssi le mode defensif est active,
 * et vaut 0 sinon.
 *
 * ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/




#ifndef __VS_DEBUG_HH
#define __VS_DEBUG_HH

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdarg.h>

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#ifndef SOTH_DEBUG_MODE
#define SOTH_DEBUG_MODE 0
#endif
#ifndef SOTH_TEMPLATE_DEBUG_MODE
#define SOTH_TEMPLATE_DEBUG_MODE 0
#endif

#define SOT_COMMON_TRACES do {  \
		    va_list arg; \
		    va_start(arg,format); \
		    vsnprintf( charbuffer,SIZE,format,arg ); \
		    va_end(arg); \
		    outputbuffer << tmpbuffer.str() << charbuffer <<std::endl; \
		} while(0)

class sotDebugTrace
{
 public:
    static const int SIZE = 512;

    std::stringstream tmpbuffer;
    std::ostream& outputbuffer;
    char charbuffer[SIZE+1];
    int traceLevel;
    int traceLevelTemplate;
    int debugPrior;

    sotDebugTrace( std::ostream& os ): outputbuffer(os),debugPrior(0) {}

    inline void trace( const int level,const char* format,...)
	{ if( level+debugPrior<=traceLevel ) SOT_COMMON_TRACES; tmpbuffer.str(""); }
    inline void trace( const char* format,...){ SOT_COMMON_TRACES;  tmpbuffer.str(""); }
    inline void trace( const int level=-1 )
	{ if( level+debugPrior<=traceLevel ) outputbuffer << tmpbuffer.str(); tmpbuffer.str("");  }

    inline void traceTemplate( const int level,const char* format,...)
	{ if( level+debugPrior<=traceLevelTemplate ) SOT_COMMON_TRACES; tmpbuffer.str(""); }
    inline void traceTemplate( const char* format,...)
	{ SOT_COMMON_TRACES; tmpbuffer.str("");  }

    inline sotDebugTrace& pre( const std::ostream& dummy ) { return *this; }
    inline sotDebugTrace& pre( const std::ostream& dummy,int level )
	{ traceLevel = level; return *this; }
/*     inline sotDebugTrace& preTemplate( const std::ostream& dummy,int level )  */
/* 	{ traceLevelTemplate = level; return *this; } */


    static const char * DEBUG_FILENAME_DEFAULT;
    static void openFile( const char * filename = DEBUG_FILENAME_DEFAULT );
    static void closeFile( const char * filename = DEBUG_FILENAME_DEFAULT );

};


extern sotDebugTrace sotDEBUGFLOW;
extern sotDebugTrace sotERRORFLOW;

#ifdef SOTH_DEBUG
#define sotPREDEBUG  "% " << __FILE__ << ": " <<__FUNCTION__  \
                              << "(#" << __LINE__ << ") :\n"
#define sotPREERROR  "\t!! "<<__FILE__ << ": " <<__FUNCTION__  \
                            << "(#" << __LINE__ << ") :"

#  define sotDEBUG(level) if( (level+sotDEBUGFLOW.debugPrior>SOTH_DEBUG_MODE)||(!sotDEBUGFLOW.outputbuffer.good()) ) ;\
    else sotDEBUGFLOW.outputbuffer<<sotPREDEBUG
#  define sotDEBUGMUTE(level) if( (level+sotDEBUGFLOW.debugPrior>SOTH_DEBUG_MODE)||(!sotDEBUGFLOW.outputbuffer.good()) ) ;\
    else sotDEBUGFLOW.outputbuffer
#  define sotERROR  if(!sotDEBUGFLOW.outputbuffer.good()); else sotERRORFLOW.outputbuffer << sotPREERROR
#  define sotDEBUGF if(!sotDEBUGFLOW.outputbuffer.good()); else sotDEBUGFLOW.pre(sotDEBUGFLOW.tmpbuffer<<sotPREDEBUG,SOTH_DEBUG_MODE).trace
#  define sotERRORF if(!sotDEBUGFLOW.outputbuffer.good()); else sotERRORFLOW.pre(sotERRORFLOW.tmpbuffer<<sotPREERROR).trace
// TEMPLATE
#  define sotTDEBUG(level) if( (level+sotDEBUGFLOW.debugPrior>SOTH_TEMPLATE_DEBUG_MODE)||(!sotDEBUGFLOW.outputbuffer.good()) ) ;\
    else sotDEBUGFLOW.outputbuffer << sotPREDEBUG
#  define sotTDEBUGF  if(!sotDEBUGFLOW.outputbuffer.good()); else sotDEBUGFLOW.pre(sotDEBUGFLOW.tmpbuffer<<sotPREDEBUG,SOTH_TEMPLATE_DEBUG_MODE).trace
inline bool sotDEBUG_ENABLE( const int & level ) { return level<=SOTH_DEBUG_MODE; }
inline bool sotTDEBUG_ENABLE( const int & level ) { return level<=SOTH_TEMPLATE_DEBUG_MODE; }

class sotDEBUGPRIORclass
{
 protected:
    int previousLevel;
    sotDEBUGPRIORclass( void ) {};
 public:
 sotDEBUGPRIORclass( unsigned int prior )
     : previousLevel(sotDEBUGFLOW.debugPrior)
	{
	    std::cout <<"Was " << sotDEBUGFLOW.debugPrior<<" set to "<<prior << std::endl;
	    sotDEBUGFLOW.debugPrior=prior;
	}
    ~sotDEBUGPRIORclass( void )
	{
	    std::cout <<"Was " << sotDEBUGFLOW.debugPrior<<" set to "<<previousLevel << std::endl;
	    sotDEBUGFLOW.debugPrior=previousLevel;
	}
};

#define sotDEBUGPRIOR(a) sotDEBUGPRIORclass sotDEBUGPRIORclass##__FUNCTION__(a);


/* -------------------------------------------------------------------------- */
#else // #ifdef SOTH_DEBUG
#define sotPREERROR  "\t!! "<<__FILE__ << ": " <<__FUNCTION__  \
                            << "(#" << __LINE__ << ") :"
#  define sotDEBUG(level) if( 1 ) ; else std::cout
#  define sotDEBUGMUTE(level) if( 1 ) ; else std::cout
#  define sotERROR sotERRORFLOW.outputbuffer << sotPREERROR
#define sotDEBUGPRIOR(a) do {} while(0)
inline void sotDEBUGF( const int level,const char* format,...) { return; }
inline void sotDEBUGF( const char* format,...) { return; }
inline void sotERRORF( const int level,const char* format,...) { return; }
inline void sotERRORF( const char* format,...) { return; }
// TEMPLATE
#  define sotTDEBUG(level) if( 1 ) ; else std::cout
inline void sotTDEBUGF( const int level,const char* format,...) { return; }
inline void sotTDEBUGF( const char* format,...) { return; }
#define sotDEBUG_ENABLE(level) false
#define sotTDEBUG_ENABLE(level) false

#endif // #ifdef SOTH_DEBUG
/* -------------------------------------------------------------------------- */

#define sotDEBUGIN(level) sotDEBUG(level) << "# In {" << std::endl
#define sotDEBUGOUT(level) sotDEBUG(level) << "# Out }" << std::endl
#define sotDEBUGINOUT(level) sotDEBUG(level) << "# In/Out { }" << std::endl

#define sotTDEBUGIN(level) sotTDEBUG(level) << "# In {" << std::endl
#define sotTDEBUGOUT(level) sotTDEBUG(level) << "# Out }" << std::endl
#define sotTDEBUGINOUT(level) sotTDEBUG(level) << "# In/Out { }" << std::endl


#endif /* #ifdef __VS_DEBUG_HH */

/*
 * Local variables:
 * c-basic-offset: 4
 * End:
 */
