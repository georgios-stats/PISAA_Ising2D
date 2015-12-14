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


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

void setseedrng(unsigned long) ;
double uniformrng(void) ;
int integerrng(int,int) ;

void get_data(char[],int*,int*,int*) ;
double cost(int*,int) ;
void initialise_from_data(int*) ;

void self_adj_grid_points(double*,int,double,double) ;
void self_adj_desired_freq(double*,int,double) ;
void self_adj_theta_update(double*,int,double*,double*,int,double*,double,double*) ;
void self_adj_theta_norm(double*,double*,double*,int,double) ;

void MH_Kpoint_int(int*,double*,int,double,int,double*,int*);

void Mutation_int_Kpoint(int*,double*, int,
						double*,double*,int,double,int,double*,
						int*) ;
void Mutation_int_Gibbs(int*,double*, int,
						double*,double*,int,
						double,int,int*) ;
void Crossover_int_Kpoint(int**,double*,
						int,int,
						double*,double*,int,double,double*,
						int*,int*) ;

void update_best_value(int *z_best, double *fz_best,
						int **x, double *fx,
						int N_dimension, int N_population ){

	int n ;
	int i ;
	int n_best ;

	n_best = 0 ;
	for ( n=1 ; n <= N_population ; n++ )
		if ( fx[n] < *fz_best) {
			n_best = n ;
			*fz_best = fx[n] ;
		}

	if (n_best > 0)
		for( i=1 ; i<=N_dimension ; i++ ) z_best[i] = x[n_best][i] ;

}

double mean_vec(double *x, int d) {

	int i ;
	double x_mean ;

	x_mean = 0.0 ;

	for (i=1; i<=d; i++) x_mean += x[i] ;

	return x_mean/d ;
}

void flags_usage(void) {
	printf("Usage:\n") ;
	printf(" -ID value \n") ;
	printf(" -Niter value \n") ;
	printf(" -Npop value \n") ;
	printf(" -Data value \n") ;
	printf(" -Gwarm value \n") ;
	printf(" -Ghigh value \n") ;
	printf(" -Gpow value \n") ;
	printf(" -Hlow value \n") ;
	printf(" -Hhigh value \n") ;
	printf(" -Hsize value \n") ;
	printf(" -Hzeta value \n") ;
	printf(" -Hconst value \n") ;
	printf(" -Twarm value \n") ;
	printf(" -Tlow value \n") ;
	printf(" -Thigh value \n") ;
	printf(" -Tlow value \n") ;
	printf(" -Tpow value \n") ;
	printf(" -flags \n") ;
	exit (8) ;
}

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))

int main(int argc, char *argv[]){

	double sumv ;
	double avev ;
	double minv ;
	double maxv ;
	double un ;

	int N_empty ;

	int lab_mcmc ;
	int i ;
	int j ;
	int pop ;
	int iter ;
	int Rep ;
	int N_sample ;
	int N_thinning ;
	int N_dimension ;
	int N_population ;
	int N_iteration ;
	int Nrow ;
	int Ncol ;

	double *freq_ref ;
	double *freq_est ;
	double freq_zeta ;
	double *theta ;
	double theta_norm_const ;
	double *grid_points ;
	int grid_size ;
	double grid_high ;
	double grid_low ;

	double gain ;
	double gain_pow ;
	int gain_warmup ;
	double gain_high ;

	double temp_saa ;
	double temp_saa_pow ;
	int temp_saa_warmup ;
	double temp_saa_high ;
	double temp_saa_low ;

	double temp_seed ;
	double temp_refine ;

	double *fx ;
	int **x ;

	double fz_best ;
	int *z_best ;

	double accpr ;

	double exp_accpr_mhseed ;
	double exp_accpr_mh0 ;
	double exp_accpr_mh1 ;
	double exp_accpr_mh2 ;
	double exp_accpr_co1 ;
	double exp_accpr_mhrefine ;

	int n_mhseed ;
	int n_mh0 ;
	int n_mh1 ;
	int n_mh2 ;
	int n_co1 ;
	int n_mhrefine ;

	int *workvec1 ;
	int *workvec2 ;
	double *workvec3 ;

	FILE *ins_hist = NULL ;
	FILE *ins_fz_best_trace = NULL ;
	FILE *ins_fz_best = NULL ;
	FILE *ins_fz_better = NULL ;

	char file_name[50] ;
	char datapath[50] ;

	/*INITIALIZE THE RNG -------------------------------------------------- */
	printf("\n INITIALIZE THE RNG \n") ;

	setseedrng( (unsigned long) time(NULL) ) ;
	for ( i=1 ; i<=10 ; i++ ) un = uniformrng() ;

	/* SET DEFAULT ALGORITHMIC SETTINGS ----------------------------------- */
	printf("\n SET ALGORITHMIC SETTINGS \n") ;

	snprintf(datapath, sizeof datapath,
			"./FPphotomicro_240x320_bw.dat") ;

	N_population = 5 ;
	N_iteration = 500000 ;
	Rep = 1 ;
	N_sample = 5000 ;
	N_thinning = N_iteration/N_sample ;

	gain_warmup = 5000 ; /*large*/
	gain_pow = 0.55 ;  /*(0.5,1] small*/
	gain_high = 1.0 ;

	temp_saa_warmup = 5000 ;
	temp_saa_pow = 0.5 ;  /*0.5*/
	temp_saa_high = 1.0 ;
	temp_saa_low = pow(10.0,-3.0) ;

	temp_seed = 100.0 ;
	temp_refine = 0.0001 ;

	freq_zeta = 0.1 ;
	grid_low = -971500.5 ;
	grid_high = -826315.5 ;
	grid_size = 200 ;
	theta_norm_const = 100.0 ;

	/*PASS EXTERNAL ALGORITHMIC ARGUMENTS -------------------------------- */

	for (i = 1; i < argc; i++) {
		if (strcmp("-ID", argv[i]) == 0)					/*..REPEAT NO*/
			Rep = atoi(argv[++i]) ;
		else if (strcmp("-Data", argv[i]) == 0) 			/*..DATA*/
			snprintf(datapath, sizeof datapath, argv[++i]) ;
		else if (strcmp("-Niter", argv[i]) == 0)			/*..COUNTERS*/
		{
			N_iteration = atoi(argv[++i]) ;
			N_thinning = N_iteration/N_sample ;
		}
		else if (strcmp("-Npop", argv[i]) == 0)
			N_population = atoi(argv[++i]) ;
		else if (strcmp("-Nsam", argv[i]) == 0)
		{
			N_sample = atoi(argv[++i]) ;
			N_thinning = N_iteration/N_sample ;
		}
		else if (strcmp("-Gwarm", argv[i]) == 0)			/*..GAIN*/
			gain_warmup = atoi(argv[++i]) ;
		else if (strcmp("-Ghigh", argv[i]) == 0)
			gain_high = atof(argv[++i]) ;
		else if (strcmp("-Gpow", argv[i]) == 0)
			gain_pow = atof(argv[++i]) ;
		else if (strcmp("-Hlow", argv[i]) == 0)				/*..HIST*/
			grid_low = atof(argv[++i]) ;
		else if (strcmp("-Hhigh", argv[i]) == 0)
			grid_high = atof(argv[++i]) ;
		else if (strcmp("-Hsize", argv[i]) == 0)
			grid_size = atoi(argv[++i]) ;
		else if (strcmp("-Hzeta", argv[i]) == 0)
			freq_zeta = atof(argv[++i]) ;
		else if (strcmp("-Hconst", argv[i]) == 0)
			theta_norm_const = atof(argv[++i]) ;
		else if (strcmp("-Twarm", argv[i]) == 0) 		/*..TEMPERATURE LADER*/
			temp_saa_warmup = atoi(argv[++i]) ;
		else if (strcmp("-Tlow", argv[i]) == 0)
			temp_saa_low = atof(argv[++i]) ;
		else if (strcmp("-Thigh", argv[i]) == 0)
			temp_saa_high = atof(argv[++i]) ;
		else if (strcmp("-Tpow", argv[i]) == 0)
			temp_saa_pow = atof(argv[++i]) ;
		else if (strcmp("-flags", argv[i]) == 0)
			flags_usage() ;
   }
	/*PRINT INPUTS OF THE ALGORITHM*/
    printf( "\n PRINT INPUTS OF THE ALGORITHM \n") ;

	printf( "Rep: \t\t %d \n", Rep) ;

	printf( "N_iteration: \t\t %d \n", N_iteration) ;
	printf( "N_population: \t\t %d \n", N_population) ;
	printf( "N_sample: \t\t %d \n", N_sample) ;
	printf( "N_thinning: \t\t %d \n", N_thinning) ;

	printf( "gain_warmup: \t\t %d \n", gain_warmup) ;
	printf( "gain_high: \t\t %f \n", gain_high) ;
	printf( "gain_pow: \t\t %f \n", gain_pow) ;

	printf( "grid_low: \t\t %f \n", grid_low) ;
	printf( "grid_high: \t\t %f \n", grid_high) ;
	printf( "grid_size: \t\t %d \n", grid_size) ;
	printf( "freq_zeta: \t\t %f \n", freq_zeta) ;
	printf( "theta_norm_const: \t %f \n", theta_norm_const) ;

	printf( "temp_saa_warmup: \t %d \n", temp_saa_warmup) ;
	printf( "temp_saa_low: \t\t %f \n", temp_saa_low) ;
	printf( "temp_saa_high: \t\t %f \n", temp_saa_high) ;
	printf( "temp_saa_pow: \t\t %f \n", temp_saa_pow) ;

	printf( "temp_seed: \t\t %f \n", temp_seed) ;
	printf( "temp_refine: \t\t %f \n", temp_refine) ;

	printf( "\n") ;

	/*GET DATA -----------------------------------------------------------*/
	printf("\n GET DATA\n") ;
	get_data(datapath,&N_dimension,&Nrow,&Ncol) ;

	/*OPEN FILES -----------------------------------------------------------*/
	printf("\n OPEN FILES\n") ;

	snprintf(file_name, sizeof file_name, "./output_files/hist-n=%d-r=%d",
			N_population, Rep);
    ins_hist = fopen( file_name , "w" ) ;

    snprintf(file_name, sizeof file_name, "./output_files/fz_best_trace-n=%d-r=%d",
    		N_population, Rep);
    ins_fz_best_trace = fopen( file_name , "w" ) ;

    snprintf(file_name, sizeof file_name, "./output_files/fz_best-n=%d-r=%d",
    		N_population, Rep);
    ins_fz_best = fopen( file_name , "w" ) ;

    snprintf(file_name, sizeof file_name, "./output_files/fz_better-n=%d-r=%d",
    		N_population, Rep);
    ins_fz_better = fopen( file_name , "w" ) ;

	/* ALLOCATE SPACE FOR THE ARRAYS -------------------------------------- */
	printf("\n ALLOCATE SPACE FOR THE ARRAYS \n") ;

	grid_points = dvector(1, grid_size) ;
	freq_ref = dvector(1, grid_size+1) ;
	theta = dvector(1, grid_size+1) ;
	freq_est = dvector(1, grid_size+1) ;

	x = imatrix(1, N_population, 1, N_dimension) ;
	fx = dvector(1, N_population) ;

	z_best = ivector(1, N_dimension) ;

	workvec1 = ivector(1, N_dimension) ;  /*working vector*/
	workvec2 = ivector(1, N_dimension) ;  /*working vector*/
	workvec3 = dvector(1, grid_size+1) ;  /*working vector*/

	/* INITIALIZE THE POPULATION -------------------------------------- */
	printf("\n INITIALIZE THE POPULATION \n") ;

	/*Get a random initialization*/

	for (pop=1 ; pop<=N_population; pop++) {
		/*for (i=1; i<=N_dimension; i++) {
			un = uniformrng() ;
			x[pop][i] = ( (un>0.5) ? 1 : 0 ) ;
		}*/
		initialise_from_data(x[pop]) ;
		fx[pop] = cost(x[pop], N_dimension) ;
	}

	/* Refine the random initialization by MRW */

	exp_accpr_mhseed = 0.0 ;  /*counter*/
	n_mhseed = 0 ; /*counter*/
	/*for (j=1 ; j<=1000 ; j++){
		temp_seed = 100.0/j ;
		for ( pop=1 ; pop<=N_population; pop++){
			MH_Kpoint_int( x[pop], &fx[pop],
							N_dimension, temp_seed,
							5, &accpr, workvec1 ) ;
			exp_accpr_mhseed += accpr ;  counter
			n_mhseed++ ;  counter
		}
	}*/
	/*exp_accpr_mhseed /= n_mhseed ;*/

	/* INITIALIZE BEST VALUES ----------------------------------------- */
	printf("\n INITIALIZE BEST VALUES \n") ;

	j = 1 ; fz_best = fx[1] ;
	for ( pop=2 ; pop <= N_population ; pop++ )
		if ( fx[pop] < fz_best ){
			j = pop ;
			fz_best = fx[j] ;
		}
	for ( i=1 ; i<=N_dimension ; i++) z_best[i] = x[j][i] ;

	/*PRINT*/

	printf( "%f \n", fz_best) ;
	for (i=1; i<=Nrow; i++) {
		for (j=1; j<=Ncol; j++) {
			if (z_best[(i-1)*Ncol+j]==1) printf( " ") ;
			if (z_best[(i-1)*Ncol+j]==0) printf( "X") ;
		}
		printf( "\n") ;
	}
	printf( "\n") ;

	/* INITIALIZE THE SELF-ADJUSTING MECHANISM ------------------------- */
	printf("\n INITIALIZE THE SELF-ADJUSTING MECHANISM \n") ;

	/* INITIALIZE THE REFERENCE FREQUENCE */
	printf("\n INITIALIZE THE REFERENCE FREQUENCE \n") ;

	self_adj_desired_freq(freq_ref, grid_size, freq_zeta) ;

	/* INITIALIZE EMPIRICAL FREQUENCE */
	printf("\n INITIALIZE EMPIRICAL FREQUENCE \n") ;

	for (i=1 ; i<=grid_size+1 ; i++) freq_est[i] = 0.0 ;

	/* INITIALIZE GRID POINTS */
	printf("\n INITIALIZE GRID POINTSE \n") ;

	self_adj_grid_points(grid_points,grid_size,grid_low,grid_high) ;

	/* INITIALIZE THETA */
	printf("\n INITIALIZE THETA \n") ;
	for (i=1 ; i<=grid_size+1 ; i++) theta[i] = 0.0 ;

	/* 	PRINT INITIAL SETTINGS OF THE SELF-ADJUSTMENT PROSEDURE */
	if( grid_size < 500 ) {
		printf("\n PRINT INITIAL SETTINGS OF "
				"THE SELF-ADJUSTMENT PROSEDURE \n") ;
		for (i=1 ; i<=grid_size ; i++)
			printf("%d \t %f \t %f \t %f \t %f \t %f \n",
					i,
					grid_points[i],
					theta[i], exp(theta[i]),
					freq_ref[i], freq_est[i] ) ;
		i = grid_size+1 ;
		printf("%d \t %f \t %f \t %f \t %f \t %f \n",
				i,
				0.0,
				theta[i], exp(theta[i]),
				freq_ref[i], freq_est[i] ) ;
	}

	/* INITIALIZE THE COUNTERS ---------------------------------------- */
	printf("\n INITIALIZE THE COUNTERS \n") ;

	/* INITIALIZE THE EXPETED ACCEPTANCE PROBABLITY COUNTERS */

	exp_accpr_mh0 = 0.0 ; n_mh0 = 0 ;
	exp_accpr_mh1 = 0.0 ; n_mh1 = 0 ;
	exp_accpr_mh2 = 0.0 ; n_mh2 = 0 ;
	exp_accpr_co1 = 0.0 ; n_co1 = 0 ;

	/* PERFORM PISAA ITERATIONS ------------------------------------------*/
	printf("\n PERFORM PISAA ITERATIONS \n") ;

	for( iter=1 ; iter<=N_iteration ; iter++ ) {

		/* UPDATE THE TEMPERATURE -------------------------------------- */

		if( iter <= temp_saa_warmup )
			temp_saa = temp_saa_high + temp_saa_low ;
		else{
			temp_saa = temp_saa_warmup/((double) iter) ;
			temp_saa = temp_saa_high*pow(temp_saa,temp_saa_pow)
							+ temp_saa_low ;
		}

		/* SAMPLING UPDATE -------------------------------------------- */

		if (N_population==1) lab_mcmc = integerrng(1, 2) ;
		else lab_mcmc = integerrng(1, 3) ;

		switch ( lab_mcmc )
		{
		case 0 :
			/*sample*/
			sumv = 0.0 ;
			for (pop=1 ; pop<=N_population ; pop++ ) {
				MH_Kpoint_int( x[pop], &fx[pop],
								N_dimension, temp_saa,
								5, &accpr, workvec1 ) ;
				sumv += accpr ;
			}
			accpr = sumv/pop ;
			exp_accpr_mh0 += accpr ;
			n_mh0++ ;
			break ;
		case 1 :
			/*sample*/
			sumv = 0.0 ;
			for (pop=1 ; pop<=N_population ; pop++ ) {
				Mutation_int_Kpoint(x[pop],&fx[pop], N_dimension,
									theta, grid_points, grid_size,
									temp_saa, 2, &accpr,
									workvec1) ;
				sumv += accpr ;
			}
			accpr = sumv/pop ;
			exp_accpr_mh1 += accpr ;
			n_mh1++ ;
			break ;
		case 2 :
			/*sample*/
			for (pop=1 ; pop<=N_population ; pop++ ) {
				Mutation_int_Gibbs(x[pop],&fx[pop], N_dimension,
									theta, grid_points, grid_size,
									temp_saa, MAX(5,MIN(25,N_population)),
									workvec1) ;
			}
			exp_accpr_mh2 += 1.0 ;
			n_mh2++ ;
			break ;
		case 3 :
			/*sample*/
			Crossover_int_Kpoint(x, fx,
									N_dimension, N_population,
									theta, grid_points, grid_size,
									temp_saa, &accpr,
									workvec1,workvec2) ;
			exp_accpr_co1 += accpr ;
			n_co1++ ;
			break ;
		default :
			break ;
		}

		/* UPDATE THE GAIN FUNCTION ------------------------------------ */

		if( iter <= gain_warmup )
			gain = gain_high ;
		else {
			gain = gain_warmup/((double) iter) ;
			gain = gain_high*pow( gain,gain_pow) ;
		}

		/* WEIGHTING UPDATE -------------------------------------------- */

		self_adj_theta_update(fx, N_population,
							theta, grid_points, grid_size,
							freq_ref, gain, workvec3) ;

		for (i=1;i<=grid_size+1;i++) freq_est[i] += workvec3[i] ;

		/* BEST VALUE UPDATE ------------------------------------------- */

		update_best_value( z_best, &fz_best, x, fx,
							N_dimension, N_population ) ;

		/* KEEP RECORDS ------------------------------------------------ */

		if ( iter > 0 && iter%N_thinning == 0  && N_sample>0)
					printf("iter=%d gain=%f temp_saa=%f min=%f \n",
							iter, gain, temp_saa, fz_best ) ;

		if (ins_fz_best_trace != NULL)
		if ( iter>0 && iter%N_thinning==0  && N_sample>0) {
			fprintf(ins_fz_best_trace,"%d %f \n", iter, fz_best) ;
			fflush(ins_fz_best_trace) ;
		}

	} /* ... ITERATIONS*/

	/*THETA-NORMALIZATION STEP ----------------------------------------- */
	printf("\n THETA-NORMALIZATION STEP \n") ;

	for ( i=1 ; i<=grid_size+1 ; i++ ) freq_est[i] /= N_iteration ;

	self_adj_theta_norm(theta,
					freq_ref,freq_est,grid_size,theta_norm_const) ;

	/* RECORD THE SELF-ADJUSTMENT MECHANISM */

	if (ins_hist != NULL){

		for ( i=1 ; i<=grid_size ; i++ )
			fprintf(ins_hist, "%d %f %f %f %f %f \n",
					i,grid_points[i],
					theta[i], exp(theta[i]),
					freq_ref[i], freq_est[i]) ;

		i = grid_size+1 ;
		fprintf(ins_hist, "%d %f %f %f %f %f \n",
				i,0.0,
				theta[i], exp(theta[i]),
				freq_ref[i], freq_est[i]) ;

		fflush(ins_hist) ;
	}

	/* RECORD BEST VALUES ---------------------------------------------- */
	printf("\n RECORD BEST VALUES \n") ;

	if (ins_fz_best != NULL){

		fprintf(ins_fz_best, "%f \n", fz_best) ;
		for (i=1; i<=Nrow; i++) {
			for (j=1; j<=Ncol; j++)
				fprintf(ins_fz_best, "%d ", z_best[(i-1)*Ncol+j]) ;
			fprintf(ins_fz_best, "\n") ;
		}
		fprintf(ins_fz_best, "\n") ;

		fprintf(ins_fz_best, "\n") ;

		fflush(ins_fz_best) ;

	}

	/* RECORD BETTER VALUES -------------------------------------------- */
	printf("\n RECORD BETTER VALUES \n") ;

	if (ins_fz_better != NULL){

		/*find*/

		maxv = fz_best ;
		for ( i=1 ; i<=N_dimension ; i++ ) workvec1[i] = z_best[i] ;
		exp_accpr_mhrefine = 0.0 ; n_mhrefine = 0 ;
		for (i=1 ; i<=10000 ; i++){
			MH_Kpoint_int(workvec1, &maxv,
							N_dimension, temp_refine,
							5, &accpr, workvec2 ) ;
			exp_accpr_mhrefine += accpr ;
			n_mhrefine++ ;
		}
		exp_accpr_mhrefine /= n_mhrefine ;

		/*print*/

		fprintf(ins_fz_better, "%f \n", maxv) ;

		for (i=1; i<=Nrow; i++) {
			for (j=1; j<=Ncol; j++)
				fprintf(ins_fz_better, "%d ", workvec1[(i-1)*Ncol+j]) ;
			fprintf(ins_fz_better, "\n") ;
		}
		fprintf(ins_fz_better, "\n") ;

		fprintf(ins_fz_better, "\n") ;

		fflush(ins_fz_better) ;

	}

	/* SUMMARY --------------------------------------------------------- */

	printf("\n SUMMARY \n") ;

	/*Result*/
	printf( "\n ...Result \n");

	printf( "fz_best : %f\n",fz_best);
	printf( "z_best  : \n");
	for (i=1; i<=Nrow; i++) {
		for (j=1; j<=Ncol; j++) {
			if (z_best[(i-1)*Ncol+j]==1) printf(" ") ;
			if (z_best[(i-1)*Ncol+j]==0) printf("X") ;
		}
		printf( "\n") ;
	}
	printf( "\n") ;

	printf( "fz_better : %f\n",maxv);
	printf( "z_better  : \n");
	for (i=1; i<=Nrow; i++) {
		for (j=1; j<=Ncol; j++) {
			if (workvec1[(i-1)*Ncol+j]==1) printf(" ") ;
			if (workvec1[(i-1)*Ncol+j]==0) printf("X") ;
		}
		printf( "\n") ;
	}
	printf( "\n") ;

	/*Acceptance ratios*/
	printf( "\n ...Acceptance ratios \n");

	exp_accpr_mh0 /= n_mh0 ;
	exp_accpr_mh1 /= n_mh1 ;
	exp_accpr_mh2 /= n_mh2 ;
	exp_accpr_co1 /= n_co1 ;

	printf( "mhs rate=%f \n",exp_accpr_mhseed) ;
	printf( "mh0 rate=%f \n",exp_accpr_mh0);
	printf( "mh1 rate=%f \n",exp_accpr_mh1);
	printf( "mh2 rate=%f \n",exp_accpr_mh2);
	printf( "co1 rate=%f \n",exp_accpr_co1);
	printf( "mhr rate=%f \n",exp_accpr_mhrefine);

	/* CLOSE FILES -------------------------------------------------------- */

	printf("\n CLOSE FILES \n") ;

	if (ins_hist != NULL) fclose( ins_hist ) ;
	if (ins_fz_best_trace != NULL) fclose( ins_fz_best_trace ) ;
	if (ins_fz_best != NULL) fclose( ins_fz_best ) ;
	if (ins_fz_better != NULL) fclose( ins_fz_better ) ;

	/* FREE SPACE --------------------------------------------------------- */

	printf("\n FREE SPACE \n") ;

	free_dvector( grid_points, 1, grid_size ) ;
	free_dvector( freq_ref, 1, grid_size+1 ) ;
	free_dvector( freq_est, 1, grid_size+1 ) ;
	free_dvector( theta, 1, grid_size+1 ) ;

	free_dvector( fx, 1, N_population ) ;
	free_imatrix( x, 1, N_population, 1, N_dimension ) ;

	free_ivector( z_best, 1, N_dimension ) ;

	free_ivector( workvec1, 1, N_dimension ) ;
	free_ivector( workvec2, 1, N_dimension ) ;
	free_dvector( workvec3, 1, grid_size+1 ) ;

	printf("\n\n\n OKAY \n") ;

	return 0 ;
}






