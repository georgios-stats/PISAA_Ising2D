/*
 * Copyrigtht 2014 Georgios Karagiannis
 *
 * This file is part of PISAA_Ising2D.
 *
 * PISAA_Ising2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 *
 * PISAA_Ising2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISAA_Ising2D.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Georgios Karagiannis 
 * Postdoctoral research associate
 * Department of Mathematics, Purdue University
 * 150 N. University Street
 * West Lafayette, IN 47907-2067, USA
 *
 * Telephone: +1 (765) 496-1007
 *
 * Email: gkaragia@purdue.edu
 *
 * Contact email: georgios.stats@gmail.com
*/


#include <math.h>

#include "RNG.h"
#include "cost_ising2D.h"
#include "MH_int_updates.h"

void MH_Kpoint_int(int *z, double *fz,
					int N_dimension, double temp,
					int K, double *accpr_pop, int *z_new) {

	int i ;
	int lab ;
	int lab2 ;
	double accpr ;
	double fy ;
	double r ;
	double un ;

    for(i=1; i<=N_dimension; i++) z_new[i] = z[i] ;

	if ( K==2 ) {
		lab = integerrng(1,N_dimension) ;
		z_new[lab] = 1-z_new[lab] ;
		lab2 = integerrng(1,N_dimension-1) ;
		if ( lab2>=lab ) i++ ;
		z_new[lab2] = 1-z_new[lab2] ;
	} else {
		for (i=1; i<=K; i++){
			lab = integerrng(1,N_dimension) ;
			z_new[lab] = 1-z_new[lab] ;
		}
	}

	fy = cost(z_new, N_dimension) ;

	r = (*fz-fy) / temp ;
	accpr = ( r>0.0 ? 1.0 : exp(r) ) ;
	un = uniformrng() ;
	if ( accpr>un ) {
		*fz = fy;
		for(i=1 ; i<=N_dimension ; i++) z[i] = z_new[i] ;
	}

	*accpr_pop = accpr ;

}

void Gibbs_int(int *z, double *fz,
		int N_dimension, double temp, int Krep, int *z_new){

	int k_old ;
	int k_new ;
	int i ;
	double Pr ;
	int lab ;
	int gg ;

	double fz_new ;
	double rat ;
	double un ;

	/* Initialize the working vector */
	for (i=1; i<=N_dimension; i++) z_new[i] = z[i] ;

	/* Perform Gibbs sweep */
	for (gg=1; gg<=Krep; gg++) {

		lab = integerrng(1, N_dimension) ;

		/* propose a new solution */

		z_new[lab] = 1 -z_new[lab] ;

		fz_new = cost(z_new, N_dimension) ;

		/* Accept */

		rat = +fz_new/temp  -*fz/temp ;

		Pr = 1.0/(1.0+exp(rat)) ;

		un = uniformrng() ;
		if ( Pr>=un ){
			*fz = fz_new;
			z[lab] = z_new[lab] ;
		} else {
			z_new[lab] = 1 -z_new[lab] ;
		}

	}

}
