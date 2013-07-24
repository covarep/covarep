
//-------------------------------------------------------------------
//   C-MEX implementation of SUMSKIPNAN - this function is part of the NaN-toolbox. 
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
// sumskipnan: sums all non-NaN values
// usage:
//	[o,count,SSQ] = sumskipnan_mex(x,DIM,flag,W);
//
// SUMSKIPNAN uses two techniques to reduce errors: 
// 1) long double (80bit) instead of 64-bit double is used internally
// 2) The Kahan Summation formula is used to reduce the error margin from N*eps to 2*eps 
//        The latter is only implemented in case of stride=1 (column vectors only, summation along 1st dimension). 
//
// Input:
// - x data array
// - DIM (optional) dimension to sum
// - flag (optional) is actually an output argument telling whether some NaN was observed
// - W (optional) weight vector to compute weighted sum (default 1)
//
// Output:
// - o (weighted) sum along dimension DIM
// - count of valid elements
// - sums of squares
//
//
//    $Id: sumskipnan_mex.cpp 9061 2011-11-11 07:47:37Z schloegl $
//    Copyright (C) 2009,2010,2011 Alois Schloegl <alois.schloegl@gmail.com>
//    This function is part of the NaN-toolbox
//    http://pub.ist.ac.at/~schloegl/matlab/NaN/
//
//-------------------------------------------------------------------




#include <math.h>
#include <stdint.h>
#include "mex.h"

inline int __sumskipnan2w__(double *data, size_t Ni, double *s, double *No, char *flag_anyISNAN, double *W);
inline int __sumskipnan3w__(double *data, size_t Ni, double *s, double *s2, double *No, char *flag_anyISNAN, double *W);
inline int __sumskipnan2wr__(double *data, size_t Ni, double *s, double *No, char *flag_anyISNAN, double *W);
inline int __sumskipnan3wr__(double *data, size_t Ni, double *s, double *s2, double *No, char *flag_anyISNAN, double *W);
inline int __sumskipnan2we__(double *data, size_t Ni, double *s, double *No, char *flag_anyISNAN, double *W);
inline int __sumskipnan3we__(double *data, size_t Ni, double *s, double *s2, double *No, char *flag_anyISNAN, double *W);
inline int __sumskipnan2wer__(double *data, size_t Ni, double *s, double *No, char *flag_anyISNAN, double *W);
inline int __sumskipnan3wer__(double *data, size_t Ni, double *s, double *s2, double *No, char *flag_anyISNAN, double *W);

//#define NO_FLAG

#ifdef tmwtypes_h
  #if (MX_API_VER<=0x07020000)
    typedef int mwSize;
  #endif 
#endif 


void mexFunction(int POutputCount,  mxArray* POutput[], int PInputCount, const mxArray *PInputs[]) 
{
    	const mwSize	*SZ;	    
    	double* 	LInput;
    	double* 	LOutputSum;
    	double* 	LOutputCount;
    	double* 	LOutputSum2;
    	long double* 	LongOutputSum = NULL;
    	long double* 	LongOutputCount = NULL;
    	long double* 	LongOutputSum2 = NULL;
    	double  	x;
    	double*		W = NULL;		// weight vector 

    	mwSize		DIM = 0; 
    	mwSize		D1, D2, D3; 	// NN; 	//  	
    	mwSize    	ND, ND2;	// number of dimensions: input, output
    	mwSize		ix0, ix1, ix2;	// index to input and output
    	mwSize    	j, l;		// running indices 
    	mwSize 		*SZ2;		// size of output 	    
	char	 	flag_isNaN = 0;

	// check for proper number of input and output arguments
	if ((PInputCount <= 0) || (PInputCount > 4))
	        mexErrMsgTxt("SUMSKIPNAN.MEX requires between 1 and 4 arguments.");
	if (POutputCount > 4)
	        mexErrMsgTxt("SUMSKIPNAN.MEX has 1 to 3 output arguments.");

	// get 1st argument
	if(mxIsDouble(PInputs[0]) && !mxIsComplex(PInputs[0]))
		LInput  = mxGetPr(PInputs[0]);
	else 	
		mexErrMsgTxt("First argument must be REAL/DOUBLE.");

    	// get 2nd argument
    	if  (PInputCount > 1) {
 	       	switch (mxGetNumberOfElements(PInputs[1])) {
		case 0: x = 0.0; 		// accept empty element
			break;
		case 1: x = (mxIsNumeric(PInputs[1]) ? mxGetScalar(PInputs[1]) : -1.0); 
			break;
		default:x = -1.0;		// invalid 
		}
		if ((x < 0) || (x > 65535) || (x != floor(x))) 
			mexErrMsgTxt("Error SUMSKIPNAN.MEX: DIM-argument must be a positive integer scalar");

		DIM = (unsigned)floor(x);	
	}

	// get size 
    	ND = mxGetNumberOfDimensions(PInputs[0]);	
    	// NN = mxGetNumberOfElements(PInputs[0]);
    	SZ = mxGetDimensions(PInputs[0]);		

	// if DIM==0 (undefined), look for first dimension with more than 1 element. 
	for (j = 0; (DIM < 1) && (j < ND); j++) 
		if (SZ[j]>1) DIM = j+1;
	
	if (DIM < 1) DIM=1;		// in case DIM is still undefined 

	ND2 = (ND>DIM ? ND : DIM);	// number of dimensions of output 

	SZ2 = (mwSize*)mxCalloc(ND2, sizeof(mwSize)); // allocate memory for output size

	for (j=0; j<ND; j++)		// copy size of input;  
		SZ2[j] = SZ[j]; 	
	for (j=ND; j<ND2; j++)		// in case DIM > ND, add extra elements 1 
		SZ2[j] = 1; 	

    	for (j=0, D1=1; j<DIM-1; D1=D1*SZ2[j++]); 	// D1 is the number of elements between two elements along dimension  DIM  
	D2 = SZ2[DIM-1];		// D2 contains the size along dimension DIM 	
    	for (j=DIM, D3=1;  j<ND; D3=D3*SZ2[j++]); 	// D3 is the number of blocks containing D1*D2 elements 

	SZ2[DIM-1] = 1;		// size of output is same as size of input but SZ(DIM)=1;

    	// get weight vector for weighted sumskipnan 
       	if  (PInputCount > 3)	{
		if (!mxGetNumberOfElements(PInputs[3])) 
			; // empty weight vector - no weighting 
		else if (mxGetNumberOfElements(PInputs[3])==D2)
			W = mxGetPr(PInputs[3]);
		else
			mexErrMsgTxt("Error SUMSKIPNAN.MEX: length of weight vector does not match size of dimension");
	}

	int ACC_LEVEL = 0;
	{
		mxArray *LEVEL = NULL;
		int s = mexCallMATLAB(1, &LEVEL, 0, NULL, "flag_accuracy_level");
		if (!s) {
			ACC_LEVEL = (int) mxGetScalar(LEVEL);
			if ((D1>1) && (ACC_LEVEL>2))
				mexWarnMsgTxt("Warning: Kahan summation not supported with stride > 1 !");
		}	
		mxDestroyArray(LEVEL);
	}
	// mexPrintf("Accuracy Level=%i\n",ACC_LEVEL);

	    // create outputs
	#define TYP mxDOUBLE_CLASS

	POutput[0] = mxCreateNumericArray(ND2, SZ2, TYP, mxREAL);
	LOutputSum = mxGetPr(POutput[0]);
	if (D1!=1 && D2>0) LongOutputSum = (long double*) mxCalloc(D1*D3,sizeof(long double));
    	if (POutputCount >= 2) {
		POutput[1] = mxCreateNumericArray(ND2, SZ2, TYP, mxREAL);
		LOutputCount = mxGetPr(POutput[1]);
		if (D1!=1 && D2>0) LongOutputCount = (long double*) mxCalloc(D1*D3,sizeof(long double));
    	}
    	if (POutputCount >= 3) {
		POutput[2] = mxCreateNumericArray(ND2, SZ2, TYP, mxREAL);
        	LOutputSum2  = mxGetPr(POutput[2]);
		if (D1!=1 && D2>0) LongOutputSum2 = (long double*) mxCalloc(D1*D3,sizeof(long double));
    	}
	mxFree(SZ2);


	if (!D1 || !D2 || !D3) // zero size array
		; 	// do nothing 
	else if (D1==1) {
	    if (ACC_LEVEL<1) {
		// double accuray, naive summation, error = N*2^-52 
		switch (POutputCount) {
		case 0: 
		case 1: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				double count;
				__sumskipnan2wr__(LInput+l*D2, D2, LOutputSum+l, &count, &flag_isNaN, W);
			}
			break;
		case 2: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				__sumskipnan2wr__(LInput+l*D2, D2, LOutputSum+l, LOutputCount+l, &flag_isNaN, W);
			}
			break;
		case 3: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				__sumskipnan3wr__(LInput+l*D2, D2, LOutputSum+l, LOutputSum2+l, LOutputCount+l, &flag_isNaN, W);
			}
			break;
		}
	    }
	    else if (ACC_LEVEL==1) {
		// extended accuray, naive summation, error = N*2^-64 
		switch (POutputCount) {
		case 0: 
		case 1: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				double count;
				__sumskipnan2w__(LInput+l*D2, D2, LOutputSum+l, &count, &flag_isNaN, W);
			}
			break;
		case 2: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				__sumskipnan2w__(LInput+l*D2, D2, LOutputSum+l, LOutputCount+l, &flag_isNaN, W);
			}
			break;
		case 3: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				__sumskipnan3w__(LInput+l*D2, D2, LOutputSum+l, LOutputSum2+l, LOutputCount+l, &flag_isNaN, W);
			}
			break;
		}
	    }
	    else if (ACC_LEVEL==3) {
		// ACC_LEVEL==3: extended accuracy and Kahan Summation, error = 2^-64
		switch (POutputCount) {
		case 0: 
		case 1: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				double count;
				__sumskipnan2we__(LInput+l*D2, D2, LOutputSum+l, &count, &flag_isNaN, W);
			}
			break;
		case 2: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				__sumskipnan2we__(LInput+l*D2, D2, LOutputSum+l, LOutputCount+l, &flag_isNaN, W);
			}
			break;
		case 3: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				__sumskipnan3we__(LInput+l*D2, D2, LOutputSum+l, LOutputSum2+l, LOutputCount+l, &flag_isNaN, W);
			}
			break;
		}
	    }
	    else if (ACC_LEVEL==2) {
		// ACC_LEVEL==2: double accuracy and Kahan Summation, error = 2^-52
		switch (POutputCount) {
		case 0: 
		case 1: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				double count;
				__sumskipnan2wer__(LInput+l*D2, D2, LOutputSum+l, &count, &flag_isNaN, W);
			}
			break;
		case 2: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				__sumskipnan2wer__(LInput+l*D2, D2, LOutputSum+l, LOutputCount+l, &flag_isNaN, W);
			}
			break;
		case 3: 
			#pragma omp parallel for schedule(dynamic)
			for (l = 0; l<D3; l++) {
				__sumskipnan3wer__(LInput+l*D2, D2, LOutputSum+l, LOutputSum2+l, LOutputCount+l, &flag_isNaN, W);
			}
			break;
		}
            }
	}
	else if (POutputCount <= 1) {
		// OUTER LOOP: along dimensions > DIM
 		for (l = 0; l<D3; l++) {
			ix0 = l*D1; 	// index for output
			ix1 = ix0*D2;	// index for input 
			for (j=0; j<D2; j++) {
				// minimize cache misses 
				ix2 =   ix0;	// index for output 
				// Inner LOOP: along dimensions < DIM
				if (W) do {
					long double x = *LInput;
        				if (!isnan(x)) {
						LongOutputSum[ix2]   += W[j]*x; 
					}
#ifndef NO_FLAG
					else 
						flag_isNaN = 1; 
#endif 
					LInput++;
					ix2++;
				} while (ix2 != (l+1)*D1);
				else do {
					long double x = *LInput;
        				if (!isnan(x)) {
						LongOutputSum[ix2]   += x; 
					}
#ifndef NO_FLAG
					else 
						flag_isNaN = 1; 
#endif 
					LInput++;
					ix2++;
				} while (ix2 != (l+1)*D1);
			}	//	end for (j=	

                        /* copy to output */
               		for (j=0; j<D1; j++) {
				LOutputSum[ix0+j] = LongOutputSum[ix0+j]; 
			}
               	}		
	}

	else if (POutputCount == 2) {
		// OUTER LOOP: along dimensions > DIM
		for (l = 0; l<D3; l++) {
			ix0 = l*D1; 
			ix1 = ix0*D2;	// index for input 
			for (j=0; j<D2; j++) {
				// minimize cache misses 
				ix2 =   ix0;	// index for output 
				// Inner LOOP: along dimensions < DIM
				if (W) do {
					long double x = *LInput;
        				if (!isnan(x)) {
						LongOutputCount[ix2] += W[j]; 
						LongOutputSum[ix2]   += W[j]*x; 
					}
#ifndef NO_FLAG
					else 
						flag_isNaN = 1; 
#endif
					LInput++;
					ix2++;
				} while (ix2 != (l+1)*D1);
				else do {
					long double x = *LInput;
        				if (!isnan(x)) {
						LongOutputCount[ix2] += 1.0; 
						LongOutputSum[ix2]   += x; 
					}
#ifndef NO_FLAG
					else 
						flag_isNaN = 1; 
#endif
					LInput++;
					ix2++;
				} while (ix2 != (l+1)*D1);
			}	//	end for (j=	

                                /* copy to output */
	               	for (j=0; j<D1; j++) {
				LOutputSum[ix0+j]   = LongOutputSum[ix0+j]; 
				LOutputCount[ix0+j] = LongOutputCount[ix0+j]; 
	       		}	// 	end else 
               	}		
	}

	else if (POutputCount == 3) {
		// OUTER LOOP: along dimensions > DIM
		for (l = 0; l<D3; l++) {
			ix0 = l*D1; 
			ix1 = ix0*D2;	// index for input 
			for (j=0; j<D2; j++) {
				// minimize cache misses 
				ix2 =   ix0;	// index for output 
				// Inner LOOP: along dimensions < DIM
				if (W) do {
					long double x = *LInput;
        				if (!isnan(x)) {
						LongOutputCount[ix2] += W[j]; 
						double t = W[j]*x;
						LongOutputSum[ix2]   += t; 
						LongOutputSum2[ix2]  += x*t; 
					}
#ifndef NO_FLAG
					else 
						flag_isNaN = 1; 
#endif
					LInput++;
					ix2++;	
				} while (ix2 != (l+1)*D1);
				else do {
					long double x = *LInput;
        				if (!isnan(x)) {
						LongOutputCount[ix2] += 1.0; 
						LongOutputSum[ix2]   += x; 
						LongOutputSum2[ix2]  += x*x; 
					}
#ifndef NO_FLAG
					else 
						flag_isNaN = 1; 
#endif
					LInput++;
					ix2++;	
				} while (ix2 != (l+1)*D1);
			}	//	end for (j=	

                       /* copy to output */
	        	for (j=0; j<D1; j++) {
				LOutputSum[ix0+j]   = LongOutputSum[ix0+j]; 
				LOutputCount[ix0+j] = LongOutputCount[ix0+j]; 
				LOutputSum2[ix0+j]  = LongOutputSum2[ix0+j]; 
			}
               	}		
	}

	if (LongOutputSum) mxFree(LongOutputSum);
	if (LongOutputCount) mxFree(LongOutputCount);
	if (LongOutputSum2) mxFree(LongOutputSum2);

#ifndef NO_FLAG
	//mexPrintf("Third argument must be not empty - otherwise status whether a NaN occured or not cannot be returned.");
	/* this is a hack, the third input argument is used to return whether a NaN occured or not. 
		this requires that the input argument is a non-empty variable
	*/	
	if  (flag_isNaN && (PInputCount > 2) && mxGetNumberOfElements(PInputs[2])) {
    		// set FLAG_NANS_OCCURED 
    		switch (mxGetClassID(PInputs[2])) {
    		case mxLOGICAL_CLASS:
    		case mxCHAR_CLASS:
    		case mxINT8_CLASS:
    		case mxUINT8_CLASS:
    			*(uint8_t*)mxGetData(PInputs[2]) = 1;
    			break; 
    		case mxDOUBLE_CLASS:
    			*(double*)mxGetData(PInputs[2]) = 1.0;
    			break; 
    		case mxSINGLE_CLASS:
    			*(float*)mxGetData(PInputs[2]) = 1.0;
    			break; 
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
    		default: 
    			mexPrintf("Type of 3rd input argument not supported.");
		}
	}
#endif
}

#define stride 1 
inline int __sumskipnan2w__(double *data, size_t Ni, double *s, double *No, char *flag_anyISNAN, double *W)
{
	long double sum=0; 
	char   flag=0; 
	// LOOP  along dimension DIM
	
	double *end = data + stride*Ni; 
	if (W) {
		// with weight vector 
		long double count = 0.0;
		do {
			long double x = *data;
        		if (!isnan(x))
			{
				count += *W; 
				sum   += *W*x;
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif	

			data++;	// stride=1
			W++;
		}
		while (data < end);
	        *No = count;
	} else {
		// w/o weight vector 
		size_t countI = 0;
		do {
			long double x = *data;
        		if (!isnan(x))
			{
				countI++; 
				sum += x; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
		}
		while (data < end);
	        *No = (double)countI;
	}	
	
#ifndef NO_FLAG
	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
#endif
	*s  = sum;

}


inline int __sumskipnan3w__(double *data, size_t Ni, double *s, double *s2, double *No, char *flag_anyISNAN, double *W)
{
	long double sum=0; 
	long double msq=0; 
	char   flag=0;
	// LOOP  along dimension DIM
	
	double *end = data + stride*Ni; 
	if (W) {
		// with weight vector 
		long double count = 0.0;
		do {
			long double x = *data;
        		if (!isnan(x)) {
				count += *W;
				long double t = *W*x; 
				sum += t; 
				msq += x*t; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
			W++;
		}
		while (data < end);
	        *No = count;
	} else {	
		// w/o weight vector 
		size_t countI = 0;
		do {
			long double x = *data;
        		if (!isnan(x)) {
				countI++; 
				sum += x; 
				msq += x*x; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
		}
		while (data < end);
	        *No = (double)countI;
	}
	
#ifndef NO_FLAG
	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
#endif
	*s  = sum;
	*s2 = msq; 
}

inline int __sumskipnan2wr__(double *data, size_t Ni, double *s, double *No, char *flag_anyISNAN, double *W)
{
	double sum=0; 
	char   flag=0; 
	// LOOP  along dimension DIM
	
	double *end = data + stride*Ni; 
	if (W) {
		// with weight vector 
		double count = 0.0;
		do {
			double x = *data;
        		if (!isnan(x))
			{
				count += *W; 
				sum   += *W*x;
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif	

			data++;	// stride=1
			W++;
		}
		while (data < end);
	        *No = count;
	} else {
		// w/o weight vector 
		size_t countI = 0;
		do {
			double x = *data;
        		if (!isnan(x))
			{
				countI++; 
				sum += x; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
		}
		while (data < end);
	        *No = (double)countI;
	}	
	
#ifndef NO_FLAG
	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
#endif
	*s  = sum;

}


inline int __sumskipnan3wr__(double *data, size_t Ni, double *s, double *s2, double *No, char *flag_anyISNAN, double *W)
{
	double sum=0; 
	double msq=0; 
	char   flag=0;
	// LOOP  along dimension DIM
	
	double *end = data + stride*Ni; 
	if (W) {
		// with weight vector 
		double count = 0.0;
		do {
			double x = *data;
        		if (!isnan(x)) {
				count += *W;
				double t = *W*x; 
				sum += t; 
				msq += x*t; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
			W++;
		}
		while (data < end);
	        *No = count;
	} else {	
		// w/o weight vector 
		size_t countI = 0;
		do {
			double x = *data;
        		if (!isnan(x)) {
				countI++; 
				sum += x; 
				msq += x*x; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
		}
		while (data < end);
	        *No = (double)countI;
	}
	
#ifndef NO_FLAG
	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
#endif
	*s  = sum;
	*s2 = msq; 
}



/***************************************
        using Kahan's summation formula [1] 
        this gives more accurate results while the computational effort within the loop is about 4x as high  
        First tests show a penalty of about 40% in terms of computational time.   

        [1] David Goldberg, 
        What Every Computer Scientist Should Know About Floating-Point Arithmetic
        ACM Computing Surveys, Vol 23, No 1, March 1991. 
 ****************************************/

inline int __sumskipnan2we__(double *data, size_t Ni, double *s, double *No, char *flag_anyISNAN, double *W)
{
	long double sum=0; 
	char   flag=0; 
	// LOOP  along dimension DIM
	
	double *end = data + stride*Ni; 
	if (W) {
		// with weight vector 
		long double count = 0.0;
		long double rc=0.0, rn=0.0;
		do {
			long double x = *data;
		        long double t,y; 
        		if (!isnan(x))
			{
				//count += *W; [1]
        			y = *W-rn;
        			t = count+y;
               			rn= (t-count)-y;
		        	count= t; 

				//sum   += *W*x; [1]
        			y = *W*x-rc;
        			t = sum+y;
               			rc= (t-sum)-y;
		        	sum= t; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif	

			data++;	// stride=1
			W++;
		}
		while (data < end);
	        *No = count;
	} else {
		// w/o weight vector 
		size_t countI = 0;
		long double rc=0.0;
		do {
			long double x = *data;
		        long double t,y; 
        		if (!isnan(x))
			{
				countI++; 
				// sum += x; [1]  
        			y = x-rc;
        			t = sum+y;
               			rc= (t-sum)-y;
		        	sum= t; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
		}
		while (data < end);
	        *No = (double)countI;
	}	
	
#ifndef NO_FLAG
	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
#endif
	*s  = sum;

}


inline int __sumskipnan3we__(double *data, size_t Ni, double *s, double *s2, double *No, char *flag_anyISNAN, double *W)
{
	long double sum=0; 
	long double msq=0; 
	char   flag=0;
	// LOOP  along dimension DIM
	
	double *end = data + stride*Ni; 
	if (W) {
		// with weight vector 
		long double count = 0.0;
        	long double rc=0.0, rn=0.0, rq=0.0;
		do {
			long double x = *data;
		        long double t,y; 
        		if (!isnan(x)) {
				//count += *W; [1]
        			y = *W-rn;
        			t = count+y;
               			rn= (t-count)-y;
		        	count= t; 

				long double w = *W*x; 
				//sum   += *W*x; [1]
        			y = *W*x-rc;
        			t = sum+y;
               			rc= (t-sum)-y;
		        	sum= t; 

				// msq += x*w; 
        			y = w*x-rq;
        			t = msq+y;
               			rq= (t-msq)-y;
		        	msq= t; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
			W++;
		}
		while (data < end);
	        *No = count;
	} else {	
		// w/o weight vector 
		size_t countI = 0;
        	long double rc=0.0, rq=0.0;
		do {
			long double x = *data;
		        long double t,y; 
        		if (!isnan(x)) {
				countI++; 
				//sum   += x; [1]
        			y = x-rc;
        			t = sum+y;
               			rc= (t-sum)-y;
		        	sum= t; 

				// msq += x*x; 
        			y = x*x-rq;
        			t = msq+y;
               			rq= (t-msq)-y;
		        	msq= t; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
		}
		while (data < end);
	        *No = (double)countI;
	}
	
#ifndef NO_FLAG
	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
#endif
	*s  = sum;
	*s2 = msq; 
}

inline int __sumskipnan2wer__(double *data, size_t Ni, double *s, double *No, char *flag_anyISNAN, double *W)
{
	double sum=0; 
	char   flag=0; 
	// LOOP  along dimension DIM
	
	double *end = data + stride*Ni; 
	if (W) {
		// with weight vector 
		double count = 0.0;
		double rc=0.0, rn=0.0;
		do {
			double x = *data;
		        double t,y; 
        		if (!isnan(x))
			{
				//count += *W; [1]
        			y = *W-rn;
        			t = count+y;
               			rn= (t-count)-y;
		        	count= t; 

				//sum   += *W*x; [1]
        			y = *W*x-rc;
        			t = sum+y;
               			rc= (t-sum)-y;
		        	sum= t; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif	

			data++;	// stride=1
			W++;
		}
		while (data < end);
	        *No = count;
	} else {
		// w/o weight vector 
		size_t countI = 0;
		double rc=0.0;
		do {
			double x = *data;
		        double t,y; 
        		if (!isnan(x))
			{
				countI++; 
				// sum += x; [1]  
        			y = x-rc;
        			t = sum+y;
               			rc= (t-sum)-y;
		        	sum= t; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
		}
		while (data < end);
	        *No = (double)countI;
	}	
	
#ifndef NO_FLAG
	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
#endif
	*s  = sum;

}


inline int __sumskipnan3wer__(double *data, size_t Ni, double *s, double *s2, double *No, char *flag_anyISNAN, double *W)
{
	double sum=0; 
	double msq=0; 
	char   flag=0;
	// LOOP  along dimension DIM
	
	double *end = data + stride*Ni; 
	if (W) {
		// with weight vector 
		double count = 0.0;
        	double rc=0.0, rn=0.0, rq=0.0;
		do {
			double x = *data;
		        double t,y; 
        		if (!isnan(x)) {
				//count += *W; [1]
        			y = *W-rn;
        			t = count+y;
               			rn= (t-count)-y;
		        	count= t; 

				double w = *W*x; 
				//sum   += *W*x; [1]
        			y = *W*x-rc;
        			t = sum+y;
               			rc= (t-sum)-y;
		        	sum= t; 

				// msq += x*w; 
        			y = w*x-rq;
        			t = msq+y;
               			rq= (t-msq)-y;
		        	msq= t; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
			W++;
		}
		while (data < end);
	        *No = count;
	} else {	
		// w/o weight vector 
		size_t countI = 0;
        	double rc=0.0, rq=0.0;
		do {
			double x = *data;
		        double t,y; 
        		if (!isnan(x)) {
				countI++; 
				//sum   += x; [1]
        			y = x-rc;
        			t = sum+y;
               			rc= (t-sum)-y;
		        	sum= t; 

				// msq += x*x; 
        			y = x*x-rq;
        			t = msq+y;
               			rq= (t-msq)-y;
		        	msq= t; 
			}
#ifndef NO_FLAG
			else 
				flag = 1; 
#endif
			data++;	// stride=1
		}
		while (data < end);
	        *No = (double)countI;
	}
	
#ifndef NO_FLAG
	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
#endif
	*s  = sum;
	*s2 = msq; 
}

