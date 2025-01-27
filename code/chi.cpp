#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "chi.h"
/*
    The implemented source of chi.h.
*/

/******************************************************************************/

double chi_square_cdf ( double x, double a )

/******************************************************************************/
/*
  Purpose:

    CHI_SQUARE_CDF evaluates the Chi squared CDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 October 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, the value of the random deviate.

    Input, double A, the parameter of the distribution, usually
    the number of degrees of freedom.

    Output, double CDF, the value of the CDF.
*/
{
  double a2;
  double b2;
  double c2;
  double cdf;
  double x2;

  x2 = 0.5 * x;
  a2 = 0.0;
  b2 = 1.0;
  c2 = 0.5 * a;

  cdf = gamma_cdf ( x2, a2, b2, c2 );

  return cdf;
}


/******************************************************************************/

double gamma_cdf ( double x, double a, double b, double c )

/******************************************************************************/
/*
  Purpose:

    GAMMA_CDF evaluates the Gamma CDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 October 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the PDF.
    A <= X

    Input, double A, B, C, the parameters of the PDF.
    0.0 < B,
    0.0 < C.

    Output, double GAMMA_CDF, the value of the CDF.
*/
{
  double cdf;
  double p2;
  double x2;

  x2 = ( x - a ) / b;
  p2 = c;

  cdf = r8_gamma_inc ( p2, x2 );

  return cdf;
}


/******************************************************************************/

double r8_gamma_inc ( double p, double x )

/******************************************************************************/
/*
  Purpose:

    R8_GAMMA_INC computes the incomplete Gamma function.

  Discussion:

    GAMMA_INC(P,  0) = 0,
    GAMMA_INC(P,+oo) = 1.

    GAMMA_INC(P,X) = Integral ( 0 <= T <= X ) T^(P-1) EXP(-T) DT / GAMMA(P).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2004

  Author:

    Original FORTRAN77 version by B L Shea.
    C version by John Burkardt.

  Reference:

    B L Shea,
    Chi-squared and Incomplete Gamma Integral,
    Algorithm AS239,
    Applied Statistics,
    Volume 37, Number 3, 1988, pages 466-473.

  Parameters:

    Input, double P, the exponent parameter.
    0.0 < P.

    Input, double X, the integral limit parameter.
    If X is less than or equal to 0, the value is returned as 0.

    Output, double R8_GAMMA_INC, the value of the function.
*/
{
  double a;
  double arg;
  double b;
  double c;
  double exp_arg_min = -88.0;
  double overflow = 1.0E+37;
  double plimit = 1000.0;
  double pn1;
  double pn2;
  double pn3;
  double pn4;
  double pn5;
  double pn6;
  double rn;
  double value;
  double tol = 1.0E-07;
  double xbig = 1.0E+08;

  value = 0.0;

  if ( p <= 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_GAMMA_INC - Fatal error!\n" );
    fprintf ( stderr, "  Parameter P <= 0.\n" );
    exit ( 1 );
  }

  if ( x <= 0.0 )
  {
    value = 0.0;
    return value;
  }
/*
  Use a normal approximation if PLIMIT < P.
*/
  if ( plimit < p )
  {
    pn1 = 3.0 * sqrt ( p ) * ( pow ( x / p, 1.0 / 3.0 )
      + 1.0 / ( 9.0 * p ) - 1.0 );
    value = normal_01_cdf ( pn1 );
    return value;
  }
/*
  Is X extremely large compared to P?
*/
  if ( xbig < x )
  {
    value = 1.0;
    return value;
  }
/*
  Use Pearson's series expansion.
  (P is not large enough to force overflow in the log of Gamma.
*/
  if ( x <= 1.0 || x < p )
  {
    arg = p * log ( x ) - x - lgamma ( p + 1.0 );
    c = 1.0;
    value = 1.0;
    a = p;

    for ( ; ; )
    {
      a = a + 1.0;
      c = c * x / a;
      value = value + c;

      if ( c <= tol )
      {
        break;
      }
    }

    arg = arg + log ( value );

    if ( exp_arg_min <= arg )
    {
      value = exp ( arg );
    }
    else
    {
      value = 0.0;
    }
  }
/*
  Use a continued fraction expansion.
*/
  else
  {
    arg = p * log ( x ) - x - lgamma ( p );
    a = 1.0 - p;
    b = a + x + 1.0;
    c = 0.0;
    pn1 = 1.0;
    pn2 = x;
    pn3 = x + 1.0;
    pn4 = x * b;
    value = pn3 / pn4;

    for ( ; ; )
    {
      a = a + 1.0;
      b = b + 2.0;
      c = c + 1.0;
      pn5 = b * pn3 - a * c * pn1;
      pn6 = b * pn4 - a * c * pn2;

      if ( 0.0 < fabs ( pn6 ) )
      {
        rn = pn5 / pn6;

        if ( fabs ( value - rn ) <= r8_min ( tol, tol * rn ) )
        {
          arg = arg + log ( value );

          if ( exp_arg_min <= arg )
          {
            value = 1.0 - exp ( arg );
          }
          else
          {
            value = 1.0;
          }

          return value;
        }
        value = rn;
      }
      pn1 = pn3;
      pn2 = pn4;
      pn3 = pn5;
      pn4 = pn6;
/*
  Rescale terms in continued fraction if terms are large.
*/
      if ( overflow <= fabs ( pn5 ) )
      {
        pn1 = pn1 / overflow;
        pn2 = pn2 / overflow;
        pn3 = pn3 / overflow;
        pn4 = pn4 / overflow;
      }
    }
  }

  return value;
}


/******************************************************************************/

double normal_01_cdf ( double x )

/******************************************************************************/
/*
  Purpose:

    NORMAL_01_CDF evaluates the Normal 01 CDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 1999

  Author:

    John Burkardt

  Reference:

    A G Adams,
    Areas Under the Normal Curve,
    Algorithm 39,
    Computer j.,
    Volume 12, pages 197-198, 1969.

  Parameters:

    Input, double X, the argument of the CDF.

    Output, double CDF, the value of the CDF.
*/
{
  double a1 = 0.398942280444;
  double a2 = 0.399903438504;
  double a3 = 5.75885480458;
  double a4 = 29.8213557808;
  double a5 = 2.62433121679;
  double a6 = 48.6959930692;
  double a7 = 5.92885724438;
  double b0 = 0.398942280385;
  double b1 = 3.8052E-08;
  double b2 = 1.00000615302;
  double b3 = 3.98064794E-04;
  double b4 = 1.98615381364;
  double b5 = 0.151679116635;
  double b6 = 5.29330324926;
  double b7 = 4.8385912808;
  double b8 = 15.1508972451;
  double b9 = 0.742380924027;
  double b10 = 30.789933034;
  double b11 = 3.99019417011;
  double cdf;
  double q;
  double y;
/*
  |X| <= 1.28.
*/
  if ( fabs ( x ) <= 1.28 )
  {
    y = 0.5 * x * x;

    q = 0.5 - fabs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5
      + a6 / ( y + a7 ) ) ) );
/*
  1.28 < |X| <= 12.7
*/
  }
  else if ( fabs ( x ) <= 12.7 )
  {
    y = 0.5 * x * x;

    q = exp ( - y ) * b0 / ( fabs ( x ) - b1
      + b2 / ( fabs ( x ) + b3
      + b4 / ( fabs ( x ) - b5
      + b6 / ( fabs ( x ) + b7
      - b8 / ( fabs ( x ) + b9
      + b10 / ( fabs ( x ) + b11 ) ) ) ) ) );
/*
  12.7 < |X|
*/
  }
  else
  {
    q = 0.0;
  }
/*
  Take account of negative X.
*/
  if ( x < 0.0 )
  {
    cdf = q;
  }
  else
  {
    cdf = 1.0 - q;
  }

  return cdf;
}


/******************************************************************************/

double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 August 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  if ( y < x )
  {
    return y;
  }
  else
  {
    return x;
  }
}
