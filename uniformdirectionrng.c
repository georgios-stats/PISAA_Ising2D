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

#include <math.h>

double normalrng(void) ;

void uniformdirectionrng(double *x, int n)
{

	double sumv ;
	int k ;

	sumv = 0 ;
	for(k = 1; k <= n; k++)
	{
	   x[k] = normalrng() ;
	   sumv += x[k]*x[k];
	}

	for(k=1; k<=n; k++) x[k] /= sqrt(sumv) ;

 }
