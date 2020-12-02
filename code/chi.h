#ifndef CHI_H_INCLUDED
#define CHI_H_INCLUDED

/*
This is a class which is utilized to calculate the CDF for chi square distribution.
*/
double chi_square_cdf ( double x, double a );
double gamma_cdf ( double x, double a, double b, double c );
double r8_gamma_inc ( double p, double x );
double normal_01_cdf ( double x );
double r8_min ( double x, double y );



#endif // CHI_H_INCLUDED
