/*
   Copyright (C)  2000    Daniel A. Atkinson  <DanAtk@aol.com>
   Copyright (C)  2004    Ivano Primi  <ivprimi@libero.it>    

   This file is part of the HPA Library.

   The HPA Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The HPA Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the HPA Library; if not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301 USA.
*/

/*
    Test:  xtan  xsin  xcos

    Uses:  atox  xprcmp  xadd  xtodbl  xprxpr

*/
#include <stdio.h>
#include "xpre.h"

int decd = 30;

int
main (void)
{
  struct xpr z, w, f, u;
  char cf[3];
  int k;
  cf[0] = 's';
  cf[1] = 'c';
  cf[2] = 't';
  for (k = 0; k < 3; ++k)
    {
      switch (cf[k])
	{
	case 't':
	  printf ("     Test of Tan Function\n");
	  break;
	case 's':
	  printf ("     Test of Sin Function\n");
	  break;
	case 'c':
	  printf ("     Test of Cos Function\n");
	  break;
	}
      z = xZero;
      w = atox ("0.2");
      u = atox ("3.01");
      for (; xprcmp (&z, &u) < 0; z = xadd (z, w, 0))
	{
	  /* compute trigonometric test function */
	  switch (cf[k])
	    {
	    case 't':
	      f = xtan (z);
	      break;
	    case 's':
	      f = xsin (z);
	      break;
	    case 'c':
	      f = xcos (z);
	      break;
	    }
	  printf (" %8.4f  ", xtodbl (z));
	  xprxpr (f, decd);
	  putchar ('\n');
	}
    }
  return 0;
}

/*  Test output

     Test of Sin Function
   0.0000  0.000000000000000000000000000000e+0
   0.2000  1.986693307950612154594126271184e-1
   0.4000  3.894183423086504916663117567957e-1
   0.6000  5.646424733950353572009454456587e-1
   0.8000  7.173560908995227616271746105814e-1
   1.0000  8.414709848078965066525023216303e-1
   1.2000  9.320390859672263496701344354948e-1
   1.4000  9.854497299884601806594745788061e-1
   1.6000  9.995736030415051643421138255462e-1
   1.8000  9.738476308781951865323731788434e-1
   2.0000  9.092974268256816953960198659117e-1
   2.2000  8.084964038195901843040369104161e-1
   2.4000  6.754631805511509265657715253413e-1
   2.6000  5.155013718214642352577269352094e-1
   2.8000  3.349881501559049195438537527124e-1
   3.0000  1.411200080598672221007448028081e-1
     Test of Cos Function
   0.0000  1.000000000000000000000000000000e+0
   0.2000  9.800665778412416311241965167482e-1
   0.4000  9.210609940028850827985267320518e-1
   0.6000  8.253356149096782972409524989554e-1
   0.8000  6.967067093471654209207499816423e-1
   1.0000  5.403023058681397174009366074430e-1
   1.2000  3.623577544766735776383733556231e-1
   1.4000  1.699671429002409386167480352036e-1
   1.6000  -2.919952230128872620577046294650e-2
   1.8000  -2.272020946930870553166743065306e-1
   2.0000  -4.161468365471423869975682295008e-1
   2.2000  -5.885011172553457085241426126549e-1
   2.4000  -7.373937155412454996088222273348e-1
   2.6000  -8.568887533689472337977021516452e-1
   2.8000  -9.422223406686581525867881173662e-1
   3.0000  -9.899924966004454572715727947313e-1
     Test of Tan Function
   0.0000  0.000000000000000000000000000000e+0
   0.2000  2.027100355086724833213582716475e-1
   0.4000  4.227932187381617619816354271653e-1
   0.6000  6.841368083416923170709254174633e-1
   0.8000  1.029638557050364012746361172820e+0
   1.0000  1.557407724654902230506974807458e+0
   1.2000  2.572151622126318935409994236033e+0
   1.4000  5.797883715482889643707720243604e+0
   1.6000  -3.423253273555741705801487543048e+1
   1.8000  -4.286261674628063525451888952280e+0
   2.0000  -2.185039863261518991643306102314e+0
   2.2000  -1.373823056768795160140036763333e+0
   2.4000  -9.160142896734105127308632475081e-1
   2.6000  -6.015966130897587227360818926913e-1
   2.8000  -3.555298316511758775773526036354e-1
   3.0000  -1.425465430742778052956354105339e-1
*/
