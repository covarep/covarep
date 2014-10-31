/* Fast linear interpolation of ordered values (see the Matlab version for doc)
%
% Copyright (c) 2012 University of Crete - Computer Science Department (UOC-CSD)
%
% License
%  This file is under the LGPL license,  you can
%  redistribute it and/or modify it under the terms of the GNU Lesser General 
%  Public License as published by the Free Software Foundation, either version 3 
%  of the License, or (at your option) any later version. This file is
%  distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
%  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
%  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
%  details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
%
% Author
%  Gilles Degottex <degottex@csd.uoc.gr>
*/


#include "mex.h"

static void interp1ordered(double *yi, double *x, double *y, double *xi, double *yid, int m, int mi)
{
/*    printf("%i:\n",__LINE__); */
    int i; /* The index of xi and yi */
    int j; /* The index of x and y */
    double g;

    i=0;
    j=0;
    while(xi[i]<x[0]){
        yi[i] = *yid;
        i++;
    }

    for(; i<mi; i++)
	{
        while(j<m-1 && x[j+1]<xi[i])
            j++;

/* printf("j=%i %i\n",j,m); */
        if(j==m-1)
            break;

        g = (xi[i]-x[j])/(x[j+1]-x[j]);
/* printf("i=%i\n",i); */
        yi[i] = y[j]*(1-g) + y[j+1]*g;
	}

	while(i<mi){
        yi[i] = *yid;
        i++;
    }
} 

 /* Entry point */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int m, mi;
    double *yi;
    double *x, *y, *xi, *yid;

   m = mxGetNumberOfElements(prhs[0]);  /* Number of input data */
   mi = mxGetNumberOfElements(prhs[2]); /* Number of data to interpolate */
   plhs[0] = mxCreateDoubleMatrix(1,mi,mxREAL); /* output vector */

   x       = mxGetPr(prhs[0]);  /* Reference abscissa */
   y       = mxGetPr(prhs[1]);  /* Reference ordinates */
   xi      = mxGetPr(prhs[2]);  /* Abscissa where to interpolate */
   yid = mxGetPr(prhs[3]);      /* Default value */
   yi = mxGetPr(plhs[0]);       /* Interpolated output values */

   interp1ordered(yi, x, y, xi, yid, m, mi);

   return;
}
