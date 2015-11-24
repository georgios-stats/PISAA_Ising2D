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
 *
 * Georgios Karagiannis © 2014  
*/

#define __COSTFUN__ 1

#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include <stdio.h>

#ifndef INFINITY
	#include <float.h>
	#define INFINITY DBL_MAX
#endif

static int *data = NULL ;
static int Ncol_data ;
static int Nrow_data ;

void get_data(char datapath[], int *N_dimension, int *Nrow, int *Ncol){

	int i ;
	int j ;
	FILE *ins;

	ins = fopen(datapath, "r") ;
	if (ins==NULL) {
		printf("No data set\n") ; abort() ;
	}

	fscanf(ins, "%d %d", &Nrow_data, &Ncol_data) ;

	data = ivector(1,Ncol_data*Nrow_data) ;

	for(i=1; i<=Nrow_data; i++)
		for(j=1; j<=Ncol_data; j++)
			fscanf(ins, "%d ", &data[(i-1)*Ncol_data+j]) ;

	*N_dimension = Ncol_data*Nrow_data ;
	*Ncol = Ncol_data ;
	*Nrow = Nrow_data ;

	fclose(ins);

	printf("\n") ;
	printf("%d %d \n",Nrow_data,Ncol_data) ;
	for(i=1; i<=Nrow_data; i++) {
		for(j=1; j<=Ncol_data; j++){
			if (data[(i-1)*Ncol_data+j]==1) printf( " ") ;
			if (data[(i-1)*Ncol_data+j]==0) printf( "X") ;
		}
		printf("\n") ;
	}
}

void initialise_from_data(int *x){
	int i ;
	int j ;
	for(i=1; i<=Nrow_data; i++)
		for(j=1; j<=Ncol_data; j++)
			x[(i-1)*Ncol_data+j] = data[(i-1)*Ncol_data+j] ;
}

#if __COSTFUN__==1

	double cost(int *x, int N_dimension){

		int i ;
		int j ;
		int ii ;
		int jj ;
		int ind ;

		double alpha_par ;
		double beta_par ;
		double gamma_par ;

		int sum0 ;
		int sum1 ;
		int sum2 ;
		double En ;

		int xx ;

		alpha_par = 1.5 ;
		beta_par = 1.5 ;
		gamma_par = 0.0 ;

		sum0 = 0; sum1 = 0 ; sum2 = 0 ;
		for (i=1; i<=Nrow_data; i++)
			for (j=1; j<=Ncol_data; j++){

				ind = (i-1)*Ncol_data+j ;

				xx = x[ ind ] ;

				sum0 += xx ;

				if (xx==data[ ind ]) sum1++ ;

				ii = i-1 ; jj = j-1 ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind = (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
				}

				ii = i-1 ; jj = j ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind = (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
				}

				ii = i-1 ; jj = j+1 ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
				}

				ii = i ; jj = j-1 ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
				}

				ii = i ; jj = j+1 ;
				if (ii>0 && jj>0) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
				}

				ii = i+1 ; jj = j-1 ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
				}

				ii = i+1 ; jj = j ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
				}

				ii = i+1 ; jj = j+1 ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
				}

			}

		En = gamma_par*( (double) sum0 )
				-alpha_par*( (double) sum1 ) -beta_par*( (double) sum2 ) ;

		return En ;

	}

#elif __COSTFUN__==2

	double cost(int *x, int N_dimension){

		int i ;
		int j ;
		int ii ;
		int jj ;
		int ind ;

		double alpha_par ;
		double beta_par ;
		double gamma_par ;

		int sum0 ;
		int sum1 ;
		int sum2 ;
		double En ;

		int xx ;

		alpha_par = 1.5 ;
		beta_par = 1.5 ;
		gamma_par = 0.0 ;

		sum0 = 0; sum1 = 0 ; sum2 = 0 ;
		for (i=1; i<=Nrow_data; i++)
			for (j=1; j<=Ncol_data; j++){

				ind = (i-1)*Ncol_data+j ;

				xx = x[ ind ] ;

				sum0 += ((xx==0)?(1):(-1)) ;

				if (xx==data[ ind ]) sum1++ ;
				else sum1-- ;

				ii = i-1 ; jj = j-1 ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
					else sum2-- ;
				}

				ii = i-1 ; jj = j ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
					else sum2-- ;
				}

				ii = i-1 ; jj = j+1 ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
					else sum2-- ;
				}

				ii = i ; jj = j-1 ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
					else sum2-- ;
				}

				ii = i ; jj = j+1 ;
				if (ii>0 && jj>0) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
					else sum2-- ;
				}

				ii = i+1 ; jj = j-1 ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
					else sum2-- ;
				}

				ii = i+1 ; jj = j ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
					else sum2-- ;
				}

				ii = i+1 ; jj = j+1 ;
				if (ii>=1 && ii<=Nrow_data && jj>=1 && jj<=Ncol_data) {
					ind =  (ii-1)*Ncol_data+jj ;
					if (xx==x[ind]) sum2++ ;
					else sum2-- ;
				}

			}

		En = gamma_par*( (double) sum0 )
				-alpha_par*( (double) sum1 ) -beta_par*( (double) sum2 ) ;

		return En ;
	}

#endif

