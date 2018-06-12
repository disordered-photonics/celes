#include "gpu/mxGPUArray.h"
#include "mex.h"
#include "wigner_D_CUDA.cuh"
#include "translation_coefficients_CUDA.cuh"

#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 256
#endif


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
@brief 	Interface to Matlab: test translation coefficients and wigner_D
=============================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	/* input arguments: 
	m, m_prime, alpha, cosBeta, gamma, kz, jTable, yTable, rResol
	*/
	
	mxInitGPU(); // initialize the MathWorks GPU API.
	
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
	
	// initialize host variables
	int const *m = (int*)mxGetData(prhs[0]);
	int const *m_prime = (int*)mxGetData(prhs[1]);
	float const *alpha = (float*)mxGetData(prhs[2]);
	float const *cosBeta = (float*)mxGetData(prhs[3]);
	float const *gamma = (float*)mxGetData(prhs[4]);
	float const *kz = (float*)mxGetData(prhs[5]);
	float const *rResol = (float*)mxGetData(prhs[8]);
	
	// initialize device variables
	float	 	*d_real_D = (float *)(mxGPUGetData(mx_real_D));
	float	 	*d_imag_D = (float *)(mxGPUGetData(mx_imag_D));
	
	float	 	*d_real_A = (float *)(mxGPUGetData(mx_real_A));
	float	 	*d_imag_A = (float *)(mxGPUGetData(mx_imag_A));
	
	float	 	*d_real_B = (float *)(mxGPUGetData(mx_real_B));
	float	 	*d_imag_B = (float *)(mxGPUGetData(mx_imag_B));
	
	float const *d_jTable = (float const *)(mxGPUGetDataReadOnly(mx_jTable));
	float const *d_yTable = (float const *)(mxGPUGetDataReadOnly(mx_yTable));
	
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
	
	mxGPUDestroyGPUArray(mx_real_D);
	mxGPUDestroyGPUArray(mx_imag_D);
	mxGPUDestroyGPUArray(mx_real_A);
	mxGPUDestroyGPUArray(mx_imag_A);
	mxGPUDestroyGPUArray(mx_real_B);
	mxGPUDestroyGPUArray(mx_imag_B);
	//cudaDeviceReset();  // necessary for profiling
}
