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
#include "Self_adjastment_prosedure.h"
#include "Mutation_int_operations.h"

/*K POINT OPERATION (INT)*/

void Mutation_int_Gibbs(int *z, double *fz, int N_dimension,
				double *theta, double *grid_points, int grid_size,
				double temp, int Krep, int *z_new){

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

		self_adj_index_search(&k_old, *fz, grid_points, grid_size) ;

		/* propose a new solution */

		z_new[lab] = 1 -z_new[lab] ;

		fz_new = cost(z_new, N_dimension) ;

		self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

		/* Accept */

		rat = theta[k_new] +fz_new/temp -theta[k_old] -*fz/temp ;

		Pr = 1.0/(1.0+exp(rat)) ;

		un = uniformrng() ;
		if ( Pr>=un ) {
			*fz = fz_new;
			z[lab] = z_new[lab] ;
		} else {
			z_new[lab] = 1 -z_new[lab] ;
		}

	}

}

/*K POINT OPERATION (INT)*/

void Mutation_int_Kpoint(int *z, double *fz, int N_dimension,
				double *theta, double *grid_points, int grid_size,
				double temp, int K, double *accpr_pop,
				int *z_new){

	/* RESTRICTION : K >= N_dimension */

	int k_old ;
	int k_new ;
	int i;
	int lab ;
	int lab2 ;

	double fz_new ;
	double rat ;
	double un ;

	self_adj_index_search(&k_old, *fz, grid_points, grid_size) ;

	/* propose a new solution */

	for (i=1; i<=N_dimension; i++) z_new[i] = z[i] ;

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

	fz_new = cost(z_new, N_dimension) ;

	self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

	rat = -theta[k_new] -fz_new/temp +theta[k_old] + *fz/temp  ;

	/* Accept reject */

	*accpr_pop = ( (rat>0.0) ? 1.0 : exp( rat ) ) ;

	un = uniformrng() ;
	if ( *accpr_pop>un ){
		*fz = fz_new;
		for ( i = 1 ; i <= N_dimension ; i++) z[i] = z_new[i] ;
	}

}







