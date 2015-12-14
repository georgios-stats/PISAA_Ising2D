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


#ifndef __SELFADJTYPE__
	#define __SELFADJTYPE__ 1
#endif

#include <math.h>

#if (__SELFADJTYPE__==1)

void self_adj_grid_points( double *grid_points, int grid_size,
						double grid_low, double grid_high ){

	double un ;
	int i ;

	un = (grid_high-grid_low)/(grid_size-1.0) ;
	for (i=1 ; i<=grid_size ; i++)
		grid_points[i] = grid_low + (i-1.0)*un ;

}

void self_adj_index_search(int *k,
					double fz, double *grid_points, int grid_size ){

	double grid_low ;
	double grid_high ;
	double grid_step ;

	grid_low = grid_points[1] ;
	grid_high = grid_points[grid_size] ;
	grid_step = (grid_high-grid_low)/(grid_size-1.0) ;

	if ( fz > grid_high ) *k = grid_size+1 ;
	else if ( fz <= grid_low ) *k = 1 ;
	else *k = (int) floor((fz-grid_low)/grid_step) +2 ;

}

#elif  (__SELFADJTYPE__==2)

void self_adj_grid_points( double *grid_points, int grid_size,
						double grid_low, double grid_high ){

	double un ;
	int i ;

	un = (grid_high-grid_low) ;

	grid_points[1] = grid_low ;

	for (i=2 ; i<=grid_size ; i++)
		grid_points[i] = grid_low + un/pow(2.0,grid_size-i) ;

}

void self_adj_index_search(int *k,
					double fz, double *grid_points, int grid_size ){

	if ( fz > grid_points[grid_size] )
		*k = grid_size+1 ;
	else
		for (*k=1; *k<=grid_size; (*k)++ )
			if (fz <= grid_points[*k]) break ;
}


#endif

void self_adj_desired_freq(double *freq_ref, int grid_size, double freq_zeta){

	int i ;
	double sumv ;

	sumv = 0.0 ;
	for (i=1 ; i<=grid_size+1 ; i++) {
		freq_ref[i] = exp( -(i-1.0)*freq_zeta ) ;
		sumv += freq_ref[i] ;
	}
	for (i=1 ; i<=grid_size+1 ; i++) freq_ref[i] /= sumv ;

}

void self_adj_theta_update(double *fx, int N_population,
						double *theta, double *grid_points, int grid_size,
						double *freq_ref, double gain, double *freq_loc){

	int n ;
	int k ;

    for ( k = 1 ; k <= grid_size+1 ; k++ ) freq_loc[k] = 0.0 ;

    for( n = 1 ; n <= N_population ; n++ ) {

    	self_adj_index_search(&k, fx[n], grid_points, grid_size ) ;

        freq_loc[k] += 1.0 ;

    }

	/* Update the histogram */

    for( k = 1 ; k <= grid_size+1 ; k++){
    	freq_loc[k] /= N_population ;
        theta[k] += gain*(freq_loc[k]-freq_ref[k]) ;
   	}

}

void self_adj_theta_norm(double *theta,
								double *freq_ref, double *freq_est,
								int grid_size, double norm_const){

	int i ;
	int N_empty ;
	double sumv ;
	double maxv ;
	double pr_e ;

	/* CORRECTION DUE TO THE EMPTY REGIONS */

	pr_e = 0.0 ;
	N_empty = 0 ;
    for ( i=1; i <= grid_size+1; i++ )
       if ( freq_est[i] <= 0.0 ){
    	   pr_e += freq_ref[i] ;
    	   N_empty++ ;
       }
    pr_e = pr_e/(grid_size+1-N_empty) ;

    for (i=1; i<=grid_size+1; i++)
    	theta[i] += log( freq_ref[i] + pr_e ) ;

    /* NORMALIZE THE BIAS PARAMETERS */

    if ( norm_const > 0.0 ) {

    	/*overflow*/
		maxv = theta[1] ;
		for	( i=1 ; i<=grid_size+1 ; i++ )
		   if ( theta[i] > maxv ) maxv = theta[i] ;

		/*stuff*/
		sumv = 0.0 ;
		for ( i=1; i<=grid_size+1; i++ ){
			theta[i] -= maxv ;
			sumv += exp(theta[i]) ;
		}

		for ( i=1; i<=grid_size+1; i++ )
			theta[i] += (-log(sumv)+log(norm_const)) ;
    }

}
