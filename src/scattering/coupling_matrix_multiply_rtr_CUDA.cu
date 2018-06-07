#include "gpu/mxGPUArray.h"
#include "mex.h"
#include "wigner_D_CUDA.cuh"
#include "translation_coefficients_CUDA.cuh"

/*=============================================================================
@brief 	Evaluate the lookup of the spherical Hankel function with cubic 
		spline interpolation

@param	p 			Spherical Hankel function order
@param	r			Radial position
@param	spjTable	Pointer to lookup table
@param	rResol		Sampling distance of radial position

@retval 			Interpolated value of spherical Hankel function
=============================================================================*/
__device__ float interpolateHankelLookup(int const p, float const r, float const *spjTable, float const rResol)
{
	float spj;
	float rPos = r/rResol;
	int rIdx = int(rPos);    					// points to table position -1, because for each p, the first entry with respect to r in the spjTable is copied 
	rPos -= rIdx; 							 	// (remainder of r/rResol) / rResol
	float rPos2 = pow(rPos,2);
	float rPos3 = pow(rPos,3);
	spj = ((-rPos3+2*rPos2-rPos) * spjTable[rIdx*(2*LMAX+1)+p]
			+ (3*rPos3-5*rPos2+2) * spjTable[(rIdx+1)*(2*LMAX+1)+p]
			+ (-3*rPos3+4*rPos2+rPos) * spjTable[(rIdx+2)*(2*LMAX+1)+p]
			+ (rPos3-rPos2) * spjTable[(rIdx+3)*(2*LMAX+1)+p])/2;
	return spj;
}


__global__ void testWignerD(int const m, int const m_prime, 
                            float const alpha, float const cosBeta, float const gamma, 
							float *real_D, float *imag_D)
{
	wignerD(m, m_prime, alpha, cosBeta, gamma, real_D, imag_D);
}


__global__ void testAxialTranslation(const float kz,
								     float const *jTable, float const *yTable,
								     const float rResol,
								     float *real_A, float *imag_A, 
								     float *real_B, float *imag_B)
{
	axialTranslationCoefficients(kz, jTable, yTable, rResol, 
	                             real_A, imag_A, real_B, imag_B);
}

/*=============================================================================
@brief 	Interface to Matlab

=============================================================================*/
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
	
	/*
	m, m_prime, alpha, cosBeta, gamma, kz, jTable, yTable, rResol
	*/
	
	
	// initialize the MathWorks GPU API.
	mxInitGPU();
	
	/*
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
	*/
	
	mxGPUArray const *mx_jTable = mxGPUCreateFromMxArray(prhs[6]);	
	mxGPUArray const *mx_yTable = mxGPUCreateFromMxArray(prhs[7]);	
	
	const mwSize dim[1] = {LMAX+1};
	mxGPUArray *mx_real_D = mxGPUCreateGPUArray((mwSize) 1,
	                                            dim,
												mxSINGLE_CLASS,
												mxREAL,
												MX_GPU_INITIALIZE_VALUES);

	mxGPUArray *mx_imag_D = mxGPUCreateGPUArray((mwSize) 1,
	                                            dim,
												mxSINGLE_CLASS,
												mxREAL,
												MX_GPU_INITIALIZE_VALUES);												

	const mwSize dim3[3] = {LMAX+1, LMAX+1, LMAX+1};
	mxGPUArray *mx_real_A = mxGPUCreateGPUArray((mwSize) 3,
	                                            dim3,
												mxSINGLE_CLASS,
												mxREAL,
												MX_GPU_INITIALIZE_VALUES);

	mxGPUArray *mx_imag_A = mxGPUCreateGPUArray((mwSize) 3,
	                                            dim3,
												mxSINGLE_CLASS,
												mxREAL,
												MX_GPU_INITIALIZE_VALUES);																								
	
	mxGPUArray *mx_real_B = mxGPUCreateGPUArray((mwSize) 3,
	                                            dim3,
												mxSINGLE_CLASS,
												mxREAL,
												MX_GPU_INITIALIZE_VALUES);

	mxGPUArray *mx_imag_B = mxGPUCreateGPUArray((mwSize) 3,
	                                            dim3,
												mxSINGLE_CLASS,
												mxREAL,
												MX_GPU_INITIALIZE_VALUES);																									
	/*
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
	*/
	
	// initialize host variables
	int const *m = (int*)mxGetData(prhs[0]);
	int const *m_prime = (int*)mxGetData(prhs[1]);
	float const *alpha = (float*)mxGetData(prhs[2]);
	float const *cosBeta = (float*)mxGetData(prhs[3]);
	float const *gamma = (float*)mxGetData(prhs[4]);
	float const *kz = (float*)mxGetData(prhs[5]);
	float const *rResol = (float*)mxGetData(prhs[8]);
	
	/*
	float const *rResol = (float*)mxGetData(prhs[9]);	// maximal polar quantum number
	int const 	threadsPerBlock = 256;
	int const blocksPerGrid = (NS[0] + threadsPerBlock - 1) / threadsPerBlock;
	*/
	
	// initialize device variables
	float	 	*d_real_D = (float *)(mxGPUGetData(mx_real_D));
	float	 	*d_imag_D = (float *)(mxGPUGetData(mx_imag_D));
	
	float	 	*d_real_A = (float *)(mxGPUGetData(mx_real_A));
	float	 	*d_imag_A = (float *)(mxGPUGetData(mx_imag_A));
	
	float	 	*d_real_B = (float *)(mxGPUGetData(mx_real_B));
	float	 	*d_imag_B = (float *)(mxGPUGetData(mx_imag_B));
	
	float const *d_jTable = (float const *)(mxGPUGetDataReadOnly(mx_jTable));
	float const *d_yTable = (float const *)(mxGPUGetDataReadOnly(mx_yTable));
	
	/*
	float const *d_real_x = (float const *)(mxGPUGetDataReadOnly(mx_real_x));
	float const *d_imag_x = (float const *)(mxGPUGetDataReadOnly(mx_imag_x));
	float const *d_real_h = (float const *)(mxGPUGetDataReadOnly(mx_real_h));
	float const *d_imag_h = (float const *)(mxGPUGetDataReadOnly(mx_imag_h));
	float const *d_PlmCoeff = (float const *)(mxGPUGetDataReadOnly(mx_PlmCoeff));
	float const *d_real_ab5 = (float const *)(mxGPUGetDataReadOnly(mx_real_ab5));
	float const *d_imag_ab5 = (float const *)(mxGPUGetDataReadOnly(mx_imag_ab5));
	float const *d_sPos = (float const *)(mxGPUGetDataReadOnly(mx_sPos));
	*/
	
	// float check_re;
	// float check_im;
	
	mexPrintf("m=");
	mexPrintf("%i\n", m[0]);
	mexPrintf("m'=");
	mexPrintf("%i\n", m_prime[0]);
	mexPrintf("alpha=");
	mexPrintf("%f\n", alpha[0]);
	mexPrintf("cos beta=");
	mexPrintf("%f\n", cosBeta[0]);
	mexPrintf("gamma=");
	mexPrintf("%f\n", gamma[0]);	
	
	testWignerD<<<1,1>>>(m[0], m_prime[0], alpha[0], cosBeta[0], gamma[0], d_real_D, d_imag_D);
	
	testAxialTranslation<<<1,1>>>(kz[0], d_jTable, d_yTable, rResol[0], 
	                              d_real_A, d_imag_A, d_real_B, d_imag_B);
	
	plhs[0] = mxGPUCreateMxArrayOnGPU(mx_real_D);
	plhs[1] = mxGPUCreateMxArrayOnGPU(mx_imag_D);
	plhs[2] = mxGPUCreateMxArrayOnGPU(mx_real_A);
	plhs[3] = mxGPUCreateMxArrayOnGPU(mx_imag_A);
	plhs[4] = mxGPUCreateMxArrayOnGPU(mx_real_B);
	plhs[5] = mxGPUCreateMxArrayOnGPU(mx_imag_B);
	
	/*
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
	*/
	
	// destroy mxgpuarrays
	/*
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
	*/
	mxGPUDestroyGPUArray(mx_real_D);
	mxGPUDestroyGPUArray(mx_imag_D);
	mxGPUDestroyGPUArray(mx_real_A);
	mxGPUDestroyGPUArray(mx_imag_A);
	mxGPUDestroyGPUArray(mx_real_B);
	mxGPUDestroyGPUArray(mx_imag_B);
	//cudaDeviceReset();  // necessary for profiling
}


