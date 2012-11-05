/*!
 * \file IDExample.h
 * \brief Header file for the abstract base class IDExample.
 * \author D. Monett D�z
 * \date September, 2006
 * 
 * This is the header of the abstract base class that defines an example.
 * Instances of this class can not be created (as it is an 'abstract class').
 * 
 * Included files:
 * 		\file adolc/adolc.h and math.h
 */
 
#ifndef IDEXAMPLE_H_
#define IDEXAMPLE_H_

#include "adolc/adolc.h"				/**< For ADOL-C functionalities. */
#include <math.h>						/**< Math. operations. */

class IDExample
{
	public:
		//
		// Constructor
		//
		IDExample( int mdyn, int mdae, double dv, char * cn )
		{
			setIndepVarTra();				// Sets nr. indep. var. trajectory, n_tra.
			setDepVarTra( mdae );			// Sets nr. dep. var. trajectory, m_tra.
			setIndepVarDyn( mdae );			// Sets nr. indep. var. dynamic, n_dyn.
			setDepVarDyn( mdyn );			// Sets nr. dep. var. dynamic, m_dyn.
			setIndepVarDAE( mdyn + mdae );	// Sets nr. indep. var. DAE, n_dae.
			setDepVarDAE( mdae );			// Sets nr. dep. var. DAE, m_dae.
			setTimeVal( dv );				// Sets value of the time.
			setName( cn );					// Sets the class name.
		}
		//
		// Destructor
		//
		//~IDExample();
		
		//
		// Pure virtual functions, to be implemented in derived clases.
		//
		virtual void tra( adouble t, adouble *ptra ) { };						/**< Trajectory. */
		virtual void dyn( adouble *px, adouble t, adouble *pd ) { };			/**< Dynamic. */
		virtual void dae( adouble *py, adouble *px, adouble t, adouble *pf ) { };/**< DAE. */
		//
		// Accessing properties
		//
		int getIndepVarTra() { return n_tra; };			/**< Returns nr. indep. var. trajectory. */
		int getDepVarTra() { return m_tra; };			/**< Returns nr. dep. var. trajectory. */
		int getIndepVarDyn() { return n_dyn; };			/**< Returns nr. indep. var. dynamic. */
		int getDepVarDyn() { return m_dyn; };			/**< Returns nr. dep. var. dynamic. */
		int getIndepVarDAE() { return n_dae; };			/**< Returns nr. indep. var. DAE. */
		int getDepVarDAE() { return m_dae; };			/**< Returns nr. dep. var. DAE. */
		double getTimeVal() { return t; };				/**< Returns value of the time. */
		char * name() { return classname; };			/**< Returns the class name. */

	private:
		int n_tra;										/**< Nr. indep. var. trajectory. */
		int m_tra;										/**< Nr. dep. var. trajectory. */
		int n_dyn;										/**< Nr. indep. var. dynamic. */
		int m_dyn;										/**< Nr. dep. var. dynamic. */
		int n_dae;										/**< Nr. indep. var. DAE. */
		int m_dae;										/**< Nr. dep. var. DAE. */
		double t;										/**< Independent var. time. */
		char * classname;								/**< Class name. */
		void setIndepVarTra() { n_tra = 1; };			/**< Sets nr. indep. var. trajectory. */
		void setDepVarTra( int iv ) { m_tra = iv; };	/**< Sets nr. dep. var. trajectory. */
		void setIndepVarDyn( int iv ) { n_dyn = iv; };	/**< Sets nr. indep. var. dynamic. */
		void setDepVarDyn( int iv ) { m_dyn = iv; };	/**< Sets nr. dep. var. dynamic. */
		void setIndepVarDAE( int iv ) { n_dae = iv; };	/**< Sets nr. indep. var. DAE. */
		void setDepVarDAE( int iv ) { m_dae = iv; };	/**< Sets nr. dep. var. DAE. */
		void setTimeVal( double dv ) { t = dv; };		/**< Sets value of the time. */
		void setName( char * cn ) { classname = cn; };	/**< Sets the class name. */
};

#endif /*IDEXAMPLE_H_*/
