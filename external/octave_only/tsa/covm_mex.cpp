/*
//-------------------------------------------------------------------
//   C-MEX implementation of COVM - this function is part of the NaN-toolbox. 
//
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, see <http://www.gnu.org/licenses/>.
//
//
// covm: in-product of matrices, NaN are skipped. 
// usage:
//	[cc,nn] = covm_mex(X,Y,flag,W);
//
// Input:
// - X:
// - Y: [optional], if empty, Y=X; 
// - flag: if not empty, it is set to 1 if some NaN was observed
// - W: weight vector to compute weighted correlation 
//
// Output:
// - CC = X' * sparse(diag(W)) * Y 	while NaN's are skipped
// - NN = real(~isnan(X)')*sparse(diag(W))*real(~isnan(Y))   count of valid (non-NaN) elements
//        computed more efficiently 
//
//    $Id: covm_mex.cpp 8223 2011-04-20 09:16:06Z schloegl $
//    Copyright (C) 2009,2010,2011 Alois Schloegl <alois.schloegl@gmail.com>
//    This function is part of the NaN-toolbox
//    http://pub.ist.ac.at/~schloegl/matlab/NaN/
//
//-------------------------------------------------------------------
*/

#ifdef __GNUC__ 
        #include <stdint.h>
#endif
#include <math.h>
#include "mex.h"

/*#define NO_FLAG*/


void mexFunction(int POutputCount,  mxArray* POutput[], int PInputCount, const mxArray *PInputs[]) 
{
    	double 		*X0=NULL, *Y0=NULL, *W=NULL;
    	double	 	*CC;
    	double 		*NN = NULL;

    	size_t		rX,cX,rY,cY;
    	size_t    	i; 
	char	 	flag_isNaN = 0;
        int             ACC_LEVEL;

	/*********** check input arguments *****************/

	// check for proper number of input and output arguments
	if ((PInputCount <= 0) || (PInputCount > 5)) {
	        mexPrintf("usage: [CC,NN] = covm_mex(X [,Y [,flag [,W [,'E']]]])\n\n");
	        mexPrintf("Do not use COVM_MEX directly, use COVM instead. \n");
/*
	        mexPrintf("\nCOVM_MEX computes the covariance matrix of real matrices and skips NaN's\n");
	        mexPrintf("\t[CC,NN] = covm_mex(...)\n\t\t computes CC=X'*Y, NN contains the number of not-NaN elements\n");
	        mexPrintf("\t\t CC./NN is the unbiased covariance matrix\n");
	        mexPrintf("\t... = covm_mex(X,Y,...)\n\t\t computes CC=X'*sparse(diag(W))*Y, number of rows of X and Y must match\n");
	        mexPrintf("\t... = covm_mex(X,[], ...)\n\t\t computes CC=X'*sparse(diag(W))*X\n");
	        mexPrintf("\t... = covm_mex(...,flag,...)\n\t\t if flag is not empty, it is set to 1 if some NaN occured in X or Y\n");
	        mexPrintf("\t... = covm_mex(...,W)\n\t\t W to compute weighted covariance, number of elements must match the number of rows of X\n");
	        mexPrintf("\t\t if isempty(W), all weights are 1\n");
	        mexPrintf("\t[CC,NN]=covm_mex(X,Y,flag,W)\n");
*/	        return;
	}
	if (POutputCount > 2)
	        mexErrMsgTxt("covm.MEX has 1 to 2 output arguments.");


	// get 1st argument
	if(mxIsDouble(PInputs[0]) && !mxIsComplex(PInputs[0]))
		X0  = mxGetPr(PInputs[0]);
	else 	
		mexErrMsgTxt("First argument must be REAL/DOUBLE.");
	rX = mxGetM(PInputs[0]);		
	cX = mxGetN(PInputs[0]);		
		
	// get 2nd argument
       	if  (PInputCount > 1)	{
		if (!mxGetNumberOfElements(PInputs[1]))
			; // Y0 = NULL; 		

		else if (mxIsDouble(PInputs[1]) && !mxIsComplex(PInputs[1]))
			Y0  = mxGetPr(PInputs[1]);
			
		else 	
			mexErrMsgTxt("Second argument must be REAL/DOUBLE.");
	}
	

    	// get weight vector for weighted sumskipnan 
       	if  (PInputCount > 3)	{
		// get 4th argument
		size_t nW = mxGetNumberOfElements(PInputs[3]);		
		if (!nW) 
			; 
		else if (nW == rX) 	
			W  = mxGetPr(PInputs[3]);
		else 	
			mexErrMsgTxt("number of elements in W must match numbers of rows in X");
	}
        
#ifdef __GNUC__
	ACC_LEVEL = 0;
	{
		mxArray *LEVEL = NULL;
		int s = mexCallMATLAB(1, &LEVEL, 0, NULL, "flag_accuracy_level");
		if (!s) {
			ACC_LEVEL = (int) mxGetScalar(LEVEL);
		}	
		mxDestroyArray(LEVEL);
	}
	// mexPrintf("Accuracy Level=%i\n",ACC_LEVEL);
#endif 
	if (Y0==NULL) {
		Y0 = X0; 
		rY = rX;
		cY = cX; 		
	}
	else {
		rY = mxGetM(PInputs[1]);		
		cY = mxGetN(PInputs[1]);		
	}
	if (rX != rY)
		mexErrMsgTxt("number of rows in X and Y do not match");

	/*********** create output arguments *****************/

	POutput[0] = mxCreateDoubleMatrix(cX, cY, mxREAL);
	CC = mxGetPr(POutput[0]);

    	if (POutputCount > 1) {
		POutput[1] = mxCreateDoubleMatrix(cX, cY, mxREAL);
		NN = mxGetPr(POutput[1]);
    	}


	/*********** compute covariance *****************/

#if 0
	/*------ version 1 --------------------- 
		this solution is slower than the alternative solution below 
		for transposed matrices, this might be faster. 
	*/	
	for (k=0; k<rX; k++) {
		double w;
		if (W) {
		w = W[k];
		for (i=0; i<cX; i++) {
			double x = X0[k+i*rX];
			if (isnan(x)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			for (j=0; j<cY; j++) {
				double y = Y0[k+j*rY];
				if (isnan(y)) {
#ifndef NO_FLAG
					flag_isNaN = 1;
#endif 
					continue;
				}
				CC[i+j*cX] += x*y*w; 
	    			if (NN != NULL) 
					NN[i+j*cX] += w; 
			}
		}
		}
		else for (i=0; i<cX; i++) {
			double x = X0[k+i*rX];
			if (isnan(x)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			for (j=0; j<cY; j++) {
				double y = Y0[k+j*rY];
				if (isnan(y)) {
#ifndef NO_FLAG
					flag_isNaN = 1;
#endif 
					continue;
				}
				CC[i+j*cX] += x*y; 
	    			if (NN != NULL) 
					NN[i+j*cX] += 1.0; 
			}
		}
	}
	
#else

   #pragma omp parallel 
   {
#ifdef __GNUC__
   if (ACC_LEVEL == 0) 
#endif
   {
	/*------ version 2 --------------------- 
	         using naive summation with double accuracy [1] 
	*/
	if ( (X0 != Y0) || (cX != cY) )
		/******** X!=Y, output is not symetric *******/	
	    if (W) /* weighted version */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		double *X = X0 + (i%cX) * rX;
		double *Y = Y0 + (i/cX) * rY;
		double cc = 0.0;
		double nw = 0.0;
		size_t k; 
		for (k=0; k<rX; k++) {
			double z = X[k]*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			cc += z*W[k];
			nw += W[k];
		}	
		CC[i] = cc; 
		if (NN != NULL) 
			NN[i] = nw; 
	    }
	    else /* no weights, all weights are 1 */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		double *X = X0 + (i%cX) * rX;
		double *Y = Y0 + (i/cX) * rY;
		double cc = 0.0;
		size_t nn = 0;
		size_t k; 
		for (k=0; k<rX; k++) {
			double z = X[k]*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			cc += z;
			nn++;
		}	
		CC[i] = cc; 
		if (NN != NULL) 
			NN[i] = (double)nn; 
	    }
	else // if (X0==Y0) && (cX==cY)
		/******** X==Y, output is symetric *******/	
	    if (W) /* weighted version */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		size_t ii = i%cX;
		size_t jj = i/cX;
	        if (ii < jj) continue;     
		double *X = X0 + ii * rX;
		double *Y = Y0 + jj * rY;
		double cc = 0.0;
		double nw = 0.0;
		size_t k; 
		for (k=0; k<rX; k++) {
			double z = X[k]*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			cc += z*W[k];
			nw += W[k];
		}
		size_t j = jj + ii*cX;
		CC[i] = cc; 
		CC[j] = cc; 
		if (NN != NULL) {
			NN[i] = nw; 
			NN[j] = nw; 
		}
	    }
	    else /* no weights, all weights are 1 */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		size_t ii = i%cX;
		size_t jj = i/cX;
	        if (ii < jj) continue;     
		double *X = X0 + ii * rX;
		double *Y = Y0 + jj * rY;
		double cc = 0.0;
		size_t nn = 0;
		size_t k; 
		for (k=0; k<rX; k++) {
			double z = X[k]*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			cc += z;
			nn++;
		}	
		size_t j = jj + ii*cX;
		CC[i] = cc; 
		CC[j] = cc; 
		if (NN != NULL) {
			NN[i] = (double)nn; 
			NN[j] = (double)nn; 
		}	
	    }

    }

#ifdef __GNUC__

    else if (ACC_LEVEL == 1) {
	/*------ version 2 --------------------- 
	         using naive summation with extended accuracy [1] 
	*/
	if ( (X0 != Y0) || (cX != cY) )
		/******** X!=Y, output is not symetric *******/	
	    if (W) /* weighted version */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		double *X = X0 + (i%cX) * rX;
		double *Y = Y0 + (i/cX) * rY;
		long double cc=0.0;
		long double nn=0.0;
		size_t k; 
		for (k=0; k<rX; k++) {
			long double z = ((long double)X[k])*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			cc += z*W[k];
			nn += W[k];
		}	
		CC[i] = cc; 
		if (NN != NULL) 
			NN[i] = nn; 
	    }
	    else /* no weights, all weights are 1 */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		double *X = X0 + (i%cX) * rX;
		double *Y = Y0 + (i/cX) * rY;
		long double cc=0.0;
		size_t nn=0;
		size_t k; 
		for (k=0; k<rX; k++) {
			long double z = ((long double)X[k])*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			cc += z;
			nn++;
		}	
		CC[i] = cc; 
		if (NN != NULL) 
			NN[i] = (double)nn; 
	    }
	else // if (X0==Y0) && (cX==cY)
		/******** X==Y, output is symetric *******/	
	    if (W) /* weighted version */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		size_t ii = i%cX;
		size_t jj = i/cX;
	        if (ii < jj) continue;     
		double *X = X0 + ii * rX;
		double *Y = Y0 + jj * rY;
		long double cc=0.0;
		long double nn=0.0;
		size_t k; 
		for (k=0; k<rX; k++) {
			long double z = ((long double)X[k])*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			cc += z*W[k];
			nn += W[k];
		}	
		size_t j = jj + ii*cX;
		CC[i] = cc; 
		CC[j] = cc; 
		if (NN != NULL) {
			NN[i] = nn; 
			NN[j] = nn; 
		}	
	    }
	    else /* no weights, all weights are 1 */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		size_t ii = i%cX;
		size_t jj = i/cX;
	        if (ii < jj) continue;     
		double *X = X0 + ii * rX;
		double *Y = Y0 + jj * rY;
		long double cc=0.0;
		size_t nn=0;
		size_t k; 
		for (k=0; k<rX; k++) {
			long double z = ((long double)X[k])*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			cc += z;
			nn++;
		}	
		size_t j = jj + ii*cX;
		CC[i] = cc; 
		CC[j] = cc; 
		if (NN != NULL) {
			NN[i] = (double)nn; 
			NN[j] = (double)nn; 
		}	
	    }

    }
    else if (ACC_LEVEL == 3) {
	/*------ version 3 --------------------- 
	        using Kahan's summation with extended (long double) accuracy [1] 
	        this gives more accurate results while the computational effort within the loop is about 4x as high  
	        However, first test show an increase in computational time of only about 25 %.   

                [1] David Goldberg, 
                What Every Computer Scientist Should Know About Floating-Point Arithmetic
                ACM Computing Surveys, Vol 23, No 1, March 1991
	*/
	if ( (X0 != Y0) || (cX != cY) )
		/******** X!=Y, output is not symetric *******/	
	    if (W) /* weighted version */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		double *X = X0 + (i%cX) * rX;
		double *Y = Y0 + (i/cX) * rY;
		long double cc=0.0;
		long double nn=0.0;
		long double rc=0.0;
		long double rn=0.0;
		size_t k; 
		for (k=0; k<rX; k++) {
		        long double t,y; 
			long double z = ((long double)X[k])*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			// cc += z*W[k]; [1]
			y = z*W[k]-rc;
			t = cc+y;
			rc= (t-cc)-y;
			cc= t; 

			// nn += W[k]; [1]
			y = z*W[k]-rn;
			t = nn+y;
			rn= (t-nn)-y;
			nn= t; 
		}	
		CC[i] = cc; 
		if (NN != NULL) 
			NN[i] = nn; 
	    }
	    else /* no weights, all weights are 1 */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		double *X = X0 + (i%cX) * rX;
		double *Y = Y0 + (i/cX) * rY;
		long double cc=0.0;
		long double rc=0.0;
		size_t nn=0;
		size_t k; 
		for (k=0; k<rX; k++) {
		        long double t,y; 
			long double z = ((long double)X[k])*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			// cc += z;  [1]
			y = z-rc;
			t = cc+y;
			rc= (t-cc)-y;
			cc= t; 

			nn++;
		}	
		CC[i] = cc; 
		if (NN != NULL) 
			NN[i] = (double)nn; 
	    }
	else // if (X0==Y0) && (cX==cY)
		/******** X==Y, output is symetric *******/	
	    if (W) /* weighted version */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		size_t ii = i%cX;
		size_t jj = i/cX;
	        if (ii < jj) continue;     
		double *X = X0 + ii * rX;
		double *Y = Y0 + jj * rY;
		long double cc=0.0;
		long double nn=0.0;
		long double rc=0.0;
		long double rn=0.0;
		size_t k; 
		for (k=0; k<rX; k++) {
		        long double t,y; 
			long double z = ((long double)X[k])*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			// cc += z*W[k];  [1]
			y = z*W[k]-rc;
			t = cc+y;
			rc= (t-cc)-y;
			cc= t; 

			// nn += W[k];  [1]
			y = z*W[k]-rn;
			t = nn+y;
			rn= (t-nn)-y;
			nn= t; 
		}	
		size_t j = jj + ii*cX;
		CC[i] = cc; 
		CC[j] = cc; 
		if (NN != NULL) {
			NN[i] = nn; 
			NN[j] = nn; 
		}	
	    }
	    else /* no weights, all weights are 1 */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		size_t ii = i%cX;
		size_t jj = i/cX;
	        if (ii < jj) continue;     
		double *X = X0 + ii * rX;
		double *Y = Y0 + jj * rY;
		long double cc=0.0;
		long double rc=0.0;
		size_t nn=0;
		size_t k; 
		for (k=0; k<rX; k++) {
		        long double t,y; 
			long double z = ((long double)X[k])*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			// cc += z; [1]
			y = z-rc;
			t = cc+y;
			rc= (t-cc)-y;
			cc= t; 

			nn++;
		}	
		size_t j = jj + ii*cX;
		CC[i] = cc; 
		CC[j] = cc; 
		if (NN != NULL) {
			NN[i] = (double)nn; 
			NN[j] = (double)nn; 
		}	
	    }
    }
    else if (ACC_LEVEL == 2) {
	/*------ version 3 --------------------- 
	        using Kahan's summation with double accuracy [1] 
	        this gives more accurate results while the computational effort within the loop is about 4x as high  
	        However, first test show an increase in computational time of only about 25 %.   

                [1] David Goldberg, 
                What Every Computer Scientist Should Know About Floating-Point Arithmetic
                ACM Computing Surveys, Vol 23, No 1, March 1991
	*/
	if ( (X0 != Y0) || (cX != cY) )
		/******** X!=Y, output is not symetric *******/	
	    if (W) /* weighted version */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		double *X = X0 + (i%cX) * rX;
		double *Y = Y0 + (i/cX) * rY;
		double cc=0.0;
		double nn=0.0;
		double rc=0.0;
		double rn=0.0;
		size_t k; 
		for (k=0; k<rX; k++) {
		        double t,y; 
			double z = X[k]*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			// cc += z*W[k]; [1]
			y = z*W[k]-rc;
			t = cc+y;
			rc= (t-cc)-y;
			cc= t; 

			// nn += W[k]; [1]
			y = z*W[k]-rn;
			t = nn+y;
			rn= (t-nn)-y;
			nn= t; 
		}	
		CC[i] = cc; 
		if (NN != NULL) 
			NN[i] = nn; 
	    }
	    else /* no weights, all weights are 1 */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		double *X = X0 + (i%cX) * rX;
		double *Y = Y0 + (i/cX) * rY;
		double cc=0.0;
		double rc=0.0;
		size_t nn=0;
		size_t k; 
		for (k=0; k<rX; k++) {
		        double t,y; 
			double z = X[k]*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			// cc += z;  [1]
			y = z-rc;
			t = cc+y;
			rc= (t-cc)-y;
			cc= t; 

			nn++;
		}	
		CC[i] = cc; 
		if (NN != NULL) 
			NN[i] = (double)nn; 
	    }
	else // if (X0==Y0) && (cX==cY)
		/******** X==Y, output is symetric *******/	
	    if (W) /* weighted version */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		size_t ii = i%cX;
		size_t jj = i/cX;
	        if (ii < jj) continue;     
		double *X = X0 + ii * rX;
		double *Y = Y0 + jj * rY;
		double cc=0.0;
		double nn=0.0;
		double rc=0.0;
		double rn=0.0;
		size_t k; 
		for (k=0; k<rX; k++) {
		        double t,y; 
			double z = X[k]*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			// cc += z*W[k];  [1]
			y = z*W[k]-rc;
			t = cc+y;
			rc= (t-cc)-y;
			cc= t; 

			// nn += W[k];  [1]
			y = z*W[k]-rn;
			t = nn+y;
			rn= (t-nn)-y;
			nn= t; 
		}	
		size_t j = jj + ii*cX;
		CC[i] = cc; 
		CC[j] = cc; 
		if (NN != NULL) {
			NN[i] = nn; 
			NN[j] = nn; 
		}	
	    }
	    else /* no weights, all weights are 1 */
            #pragma omp for schedule(dynamic) nowait
	    for (i = 0; i < cX * cY; i++)
	    {
		size_t ii = i%cX;
		size_t jj = i/cX;
	        if (ii < jj) continue;     
		double *X = X0 + ii * rX;
		double *Y = Y0 + jj * rY;
		double cc=0.0;
		double rc=0.0;
		size_t nn=0;
		size_t k; 
		for (k=0; k<rX; k++) {
		        double t,y; 
			double z = X[k]*Y[k];
			if (isnan(z)) {
#ifndef NO_FLAG
				flag_isNaN = 1;
#endif 
				continue;
			}
			// cc += z; [1]
			y = z-rc;
			t = cc+y;
			rc= (t-cc)-y;
			cc= t; 

			nn++;
		}	
		size_t j = jj + ii*cX;
		CC[i] = cc; 
		CC[j] = cc; 
		if (NN != NULL) {
			NN[i] = (double)nn; 
			NN[j] = (double)nn; 
		}	
	    }
    }
#endif
   } // end pragma omg parallel 
   
   
#ifndef NO_FLAG
	//mexPrintf("Third argument must be not empty - otherwise status whether a NaN occured or not cannot be returned.");
	/* this is a hack, the third input argument is used to return whether a NaN occured or not. 
		this requires that the input argument is a non-empty variable
	*/	
	if  (flag_isNaN && (PInputCount > 2) && mxGetNumberOfElements(PInputs[2])) {
    		// set FLAG_NANS_OCCURED 
    		switch (mxGetClassID(PInputs[2])) {
    		case mxDOUBLE_CLASS:
    			*(double*)mxGetData(PInputs[2]) = 1.0;
    			break; 
    		case mxSINGLE_CLASS:
    			*(float*)mxGetData(PInputs[2]) = 1.0;
    			break; 
    		case mxLOGICAL_CLASS:
    		case mxCHAR_CLASS:
    		case mxINT8_CLASS:
    		case mxUINT8_CLASS:
    			*(char*)mxGetData(PInputs[2]) = 1;
    			break; 
#ifdef __GNUC__ 
    		case mxINT16_CLASS:
    		case mxUINT16_CLASS:
    			*(uint16_t*)mxGetData(PInputs[2]) = 1;
    			break; 
    		case mxINT32_CLASS:
    		case mxUINT32_CLASS:
    			*(uint32_t*)mxGetData(PInputs[2])= 1;
    			break; 
    		case mxINT64_CLASS:
    		case mxUINT64_CLASS:
    			*(uint64_t*)mxGetData(PInputs[2]) = 1;
    			break; 
    		case mxFUNCTION_CLASS:
    		case mxUNKNOWN_CLASS:
    		case mxCELL_CLASS:
    		case mxSTRUCT_CLASS:
#endif
    		default: 
    			mexPrintf("Type of 3rd input argument cannot be used to return status of NaN occurence.");
		}
	}
#endif
#endif
}

