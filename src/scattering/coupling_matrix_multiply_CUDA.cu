/**
 * Copyright (c) 2017, Amos Egel (KIT), Lorenzo Pattelli (LENS)
 *                     Giacomo Mazzamuto (LENS)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 *  * Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "matrix.h"
#include <math.h>
#include "cuda_profiler_api.h"

__device__ float assocLegendreFunction(int const l, int const m, float const ct, float const st, float const *plm_coeffs)
{
	float Plm = 0.0f;
	int jj=0;
	for (int lambda=l-m; lambda>=0; lambda-=2)
	{
		Plm = Plm + pow(st,m) * pow(ct,lambda) * plm_coeffs[jj*(2*LMAX+1)*(2*LMAX+1)+m*(2*LMAX+1)+l];
		jj++;
	}
	return Plm;
}

__device__ float sphericalBesselLookup(int const p, float const r, float const *spjTable, float const rResol)
{
	float spj = 0.0f;
	float rPos = r/rResol;
	int rIdx = int(rPos);    						// points to table position -1, because for each p, the first entry with respect to r in the spjTable is copied 
	rPos -= rIdx; 							 	// (remainder of r/rResol) / rResol
	float rPos2 = pow(rPos,2);
	float rPos3 = pow(rPos,3);
	spj = ((-rPos3+2*rPos2-rPos) * spjTable[rIdx*(2*LMAX+1)+p]
			+ (3*rPos3-5*rPos2+2) * spjTable[(rIdx+1)*(2*LMAX+1)+p]
			+ (-3*rPos3+4*rPos2+rPos) * spjTable[(rIdx+2)*(2*LMAX+1)+p]
			+ (rPos3-rPos2) * spjTable[(rIdx+3)*(2*LMAX+1)+p])/2;
	return spj;
}


__global__ void translationMatrixProduct(int const s2, int const NS, float const *sPosArray,
										float const *sphericalBesselTable, float const *sphericalNeumannTable, float rResol,
										float const *plm_coeffs, float const *re_abTable, float const *im_abTable,
										float const *re_x, float const *im_x, float *re_Wx, float *im_Wx)
{
  int const s1 = blockDim.x * blockIdx.x + threadIdx.x + 1; // receiving sphere number (1...NS)
  float x21, y21, z21;
	float r, cosTheta, sinTheta, phi;
	float re_h[2*LMAX+1];
	float im_h[2*LMAX+1];
	float Ppdm[(2*LMAX+1)*(2*LMAX+2)/2];
	float cosmphi[4*LMAX+1];
	float sinmphi[4*LMAX+1];
	int n1, n2, deltam;
	float re_xTmp, im_xTmp;
	int loopCounter = 0;
	int WxIdx, xIdx, abIdx;
	float re_incr, im_incr;
	float re_abP, im_abP, re_abPh, im_abPh, re_abPheimp, im_abPheimp;
	
  if ((s1!=s2)&&s1<=NS)
	{
		// relative position
		
		x21 = sPosArray[3*(s1-1)]-sPosArray[3*(s2-1)];
		y21 = sPosArray[3*(s1-1)+1]-sPosArray[3*(s2-1)+1];
		z21 = sPosArray[3*(s1-1)+2]-sPosArray[3*(s2-1)+2];
		
		r = sqrt(x21*x21+y21*y21+z21*z21);
		cosTheta = z21/r;
		sinTheta = sqrt(1-cosTheta*cosTheta);
		phi = atan2(y21,x21);
		
		for (int p=0; p<=2*LMAX; p++)	// precompute spherical Hankel functions and Legendre functions
		{
			re_h[p] = sphericalBesselLookup(p,r,sphericalBesselTable,rResol);
			im_h[p] = sphericalBesselLookup(p,r,sphericalNeumannTable,rResol);
			for (int absdm=0; absdm<=p; absdm++)
			{
				Ppdm[p*(p+1)/2+absdm] = assocLegendreFunction(p,absdm,cosTheta,sinTheta,plm_coeffs);
			}
		}
		
		for (int dm=-2*LMAX; dm<=2*LMAX; dm++) // precompute exp(i(m-m')phi)
		{
			cosmphi[dm+2*LMAX] = cosf(dm*phi);
			sinmphi[dm+2*LMAX] = sinf(dm*phi);
		}
		
		for (int tau1=1; tau1<=2; tau1++) // evaluate matrix vector product
		{
			int temp1 = (tau1-1)*LMAX*(LMAX+2);
			for (int l1=1; l1<=LMAX; l1++)
			{
				int coeff1 = temp1+(l1-1)*(l1+1)+l1+1;
				for (int m1=-l1; m1<=l1; m1++)
				{
					n1 = coeff1+m1;
					WxIdx = (n1-1)*NS+s1-1;
					re_incr = 0.0f;
					im_incr = 0.0f;
					
					for (int tau2=1; tau2<=2; tau2++)
					{
						int temp2 = (tau2-1)*LMAX*(LMAX+2);
						for (int l2=1; l2<=LMAX; l2++)
						{
							int coeff2 = temp2+(l2-1)*(l2+1)+l2+1;
							for (int m2=-l2; m2<=l2; m2++)
							{
								n2 = coeff2+m2;
								xIdx = (n2-1)*NS+s2-1;
								re_xTmp = re_x[xIdx];
								im_xTmp = im_x[xIdx];
								deltam=m2-m1;
								abIdx=deltam+2*LMAX;
								for (int p=max(abs(deltam),abs(l1-l2)+abs(tau1-tau2)); p<=l1+l2; p++)
								{
									re_abP = re_abTable[loopCounter]*Ppdm[p*(p+1)/2+abs(deltam)];
									im_abP = im_abTable[loopCounter]*Ppdm[p*(p+1)/2+abs(deltam)];
									
									re_abPh = re_abP*re_h[p] - im_abP*im_h[p];
									im_abPh = re_abP*im_h[p] + im_abP*re_h[p];
									
									re_abPheimp = re_abPh*cosmphi[abIdx] - im_abPh*sinmphi[abIdx];
									im_abPheimp = re_abPh*sinmphi[abIdx] + im_abPh*cosmphi[abIdx];
									
									re_incr += re_abPheimp*re_xTmp - im_abPheimp*im_xTmp;
									im_incr += re_abPheimp*im_xTmp + im_abPheimp*re_xTmp;
									
									loopCounter++;
								}	
							}
						}
					}//tau2
					re_Wx[WxIdx] += re_incr;
					im_Wx[WxIdx] += im_incr;
				} 
			}  
		}//tau1
	}  
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	/* input arguments: 
	real_x
	imag_x 
	real_h			
	imag_h	
	Plm_coeffs
	real_ab5
	imag_ab5
	spos
	NS 			
	rResol */

	// initialize the MathWorks GPU API.
	mxInitGPU();
		
	// check number of arguments:
	if (nrhs!=10) {mexErrMsgTxt("wrong number of input arguments");}
	if (nlhs!=2) {mexErrMsgTxt("wrong number of output arguments");}

	// check for GPUArrays
	if (!(mxIsGPUArray(prhs[0]))) {mexErrMsgTxt("real_x is not a gpuArray");}
	if (!(mxIsGPUArray(prhs[1]))) {mexErrMsgTxt("imag_x is not a gpuArray");}
	if (!(mxIsGPUArray(prhs[2]))) {mexErrMsgTxt("real_h is not a gpuArray");}
	if (!(mxIsGPUArray(prhs[3]))) {mexErrMsgTxt("imag_h is not a gpuArray");}
	if (!(mxIsGPUArray(prhs[4]))) {mexErrMsgTxt("Plm_coeffs is not a gpuArray");}
	if (!(mxIsGPUArray(prhs[5]))) {mexErrMsgTxt("real_ab5 is not a gpuArray");}
	if (!(mxIsGPUArray(prhs[6]))) {mexErrMsgTxt("imag_ab5 is not a gpuArray");}
	if (!(mxIsGPUArray(prhs[7]))) {mexErrMsgTxt("spos is not a gpuArray");}
	
	// initialize mxGPUArrays
	mxGPUArray const *mx_real_x = mxGPUCreateFromMxArray(prhs[0]);	
	mxGPUArray const *mx_imag_x = mxGPUCreateFromMxArray(prhs[1]);	
	mxGPUArray const *mx_real_h = mxGPUCreateFromMxArray(prhs[2]);	
	mxGPUArray const *mx_imag_h = mxGPUCreateFromMxArray(prhs[3]);	
	mxGPUArray const *mx_PlmCoeff = mxGPUCreateFromMxArray(prhs[4]);	
	mxGPUArray const *mx_real_ab5 = mxGPUCreateFromMxArray(prhs[5]);
	mxGPUArray const *mx_imag_ab5 = mxGPUCreateFromMxArray(prhs[6]);
	mxGPUArray const *mx_sPos = mxGPUCreateFromMxArray(prhs[7]);
	// mxGPUArray *mx_real_Wx = mxGPUCopyGPUArray(mx_real_x);  // writable copy of real_x ... to be overwritten in kernel
	// mxGPUArray *mx_imag_Wx = mxGPUCopyGPUArray(mx_imag_x);

	mxGPUArray *mx_real_Wx = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(mx_real_x),
                            mxGPUGetDimensions(mx_real_x),
                            mxSINGLE_CLASS,mxREAL,MX_GPU_INITIALIZE_VALUES);
																
	mxGPUArray *mx_imag_Wx = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(mx_real_x),
                            mxGPUGetDimensions(mx_real_x),
                            mxSINGLE_CLASS,mxREAL,MX_GPU_INITIALIZE_VALUES);
	
	// check data types
	if (mxGPUGetClassID(mx_real_x) != mxSINGLE_CLASS) {mexErrMsgTxt("real_x is not single");}
	if (mxGPUGetClassID(mx_imag_x) != mxSINGLE_CLASS) {mexErrMsgTxt("imag_x is not single");}
	if (mxGPUGetClassID(mx_real_h) != mxSINGLE_CLASS) {mexErrMsgTxt("real_h is not single");}
	if (mxGPUGetClassID(mx_imag_h) != mxSINGLE_CLASS) {mexErrMsgTxt("imag_h is not single");}
	if (mxGPUGetClassID(mx_PlmCoeff) != mxSINGLE_CLASS) {mexErrMsgTxt("PlmCoeff is not single");}
	if (mxGPUGetClassID(mx_real_ab5) != mxSINGLE_CLASS) {mexErrMsgTxt("real_ab5 is not single");}
	if (mxGPUGetClassID(mx_imag_ab5) != mxSINGLE_CLASS) {mexErrMsgTxt("imag_ab5 is not single");}
	if (mxGPUGetClassID(mx_sPos) != mxSINGLE_CLASS) {mexErrMsgTxt("sPos is not single");}
	if (mxGetClassID(prhs[8]) != mxINT32_CLASS) {mexErrMsgTxt("NS is not int32");}
	if (mxGetClassID(prhs[9]) != mxSINGLE_CLASS) {mexErrMsgTxt("rResol is not single");}
	
	// initialize host variables
	int const *NS = (int*)mxGetData(prhs[8]);	// total number of spheres
	float const *rResol = (float*)mxGetData(prhs[9]);	// maximal polar quantum number
	int const 	threadsPerBlock = 256;
	int const blocksPerGrid = (NS[0] + threadsPerBlock - 1) / threadsPerBlock;

	// initialize device variables
	float	 	*d_real_Wx = (float *)(mxGPUGetData(mx_real_Wx));
	float	 	*d_imag_Wx = (float *)(mxGPUGetData(mx_imag_Wx));
	float const *d_real_x = (float const *)(mxGPUGetDataReadOnly(mx_real_x));
	float const *d_imag_x = (float const *)(mxGPUGetDataReadOnly(mx_imag_x));
	float const *d_real_h = (float const *)(mxGPUGetDataReadOnly(mx_real_h));
	float const *d_imag_h = (float const *)(mxGPUGetDataReadOnly(mx_imag_h));
	float const *d_PlmCoeff = (float const *)(mxGPUGetDataReadOnly(mx_PlmCoeff));
	float const *d_real_ab5 = (float const *)(mxGPUGetDataReadOnly(mx_real_ab5));
	float const *d_imag_ab5 = (float const *)(mxGPUGetDataReadOnly(mx_imag_ab5));
	float const *d_sPos = (float const *)(mxGPUGetDataReadOnly(mx_sPos));

	
	// float check_re;
	// float check_im;
	
	// start computation
	for (int s2=1; s2<=NS[0]; s2++)
	{
		translationMatrixProduct<<< blocksPerGrid,threadsPerBlock >>> (s2,NS[0],d_sPos,
																		d_real_h, d_imag_h, rResol[0],
																		d_PlmCoeff, d_real_ab5, d_imag_ab5,
																		d_real_x, d_imag_x, d_real_Wx, d_imag_Wx);
		// cudaMemcpy(&check_re,d_real_Wx,sizeof(check_re),cudaMemcpyDeviceToHost);
		// cudaMemcpy(&check_im,d_imag_Wx,sizeof(check_im),cudaMemcpyDeviceToHost);
		// mexPrintf("%f %f\n",check_re,check_im);
																																
	}

	// wrap the result up as a MATLAB gpuArray for return
	plhs[0] = mxGPUCreateMxArrayOnGPU(mx_real_Wx);
	plhs[1] = mxGPUCreateMxArrayOnGPU(mx_imag_Wx);

	// destroy mxgpuarrays
	mxGPUDestroyGPUArray(mx_real_x);
	mxGPUDestroyGPUArray(mx_imag_x);
	mxGPUDestroyGPUArray(mx_real_h);
	mxGPUDestroyGPUArray(mx_imag_h);
	mxGPUDestroyGPUArray(mx_PlmCoeff);
	mxGPUDestroyGPUArray(mx_real_ab5);
	mxGPUDestroyGPUArray(mx_imag_ab5);
	mxGPUDestroyGPUArray(mx_sPos);
	mxGPUDestroyGPUArray(mx_real_Wx);
	mxGPUDestroyGPUArray(mx_imag_Wx);
		
	//cudaDeviceReset();  // necessary for profiling
}
