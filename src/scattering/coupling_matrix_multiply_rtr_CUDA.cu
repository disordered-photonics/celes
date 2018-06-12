#include "gpu/mxGPUArray.h"
#include "mex.h"
#include "kernel_rtr.cu"

#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 128
#endif


/*=============================================================================
@brief 	Interface to Matlab: coupling matrix multiply W*x
=============================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	/* input arguments: 
	real_x
	imag_x 
	real_h			
	imag_h	
	spos
	NS 			
	rResol */
	
	// initialize the MathWorks GPU API.
	mxInitGPU();
	
	// check number of arguments:
	if (nrhs!=7) {mexErrMsgTxt("wrong number of input arguments");}
	if (nlhs!=2) {mexErrMsgTxt("wrong number of output arguments");}
	
	// check for GPUArrays
	if (!(mxIsGPUArray(prhs[0]))) {mexErrMsgTxt("real_x is not a gpuArray");}
	if (!(mxIsGPUArray(prhs[1]))) {mexErrMsgTxt("imag_x is not a gpuArray");}
	if (!(mxIsGPUArray(prhs[2]))) {mexErrMsgTxt("real_h is not a gpuArray");}
	if (!(mxIsGPUArray(prhs[3]))) {mexErrMsgTxt("imag_h is not a gpuArray");}
	if (!(mxIsGPUArray(prhs[4]))) {mexErrMsgTxt("spos is not a gpuArray");}
	
	// initialize mxGPUArrays
	mxGPUArray const *mx_real_x = mxGPUCreateFromMxArray(prhs[0]);	
	mxGPUArray const *mx_imag_x = mxGPUCreateFromMxArray(prhs[1]);	
	mxGPUArray const *mx_real_h = mxGPUCreateFromMxArray(prhs[2]);	
	mxGPUArray const *mx_imag_h = mxGPUCreateFromMxArray(prhs[3]);	
	mxGPUArray const *mx_sPos = mxGPUCreateFromMxArray(prhs[4]);
	
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
	if (mxGPUGetClassID(mx_sPos) != mxSINGLE_CLASS) {mexErrMsgTxt("sPos is not single");}
	if (mxGetClassID(prhs[5]) != mxINT32_CLASS) {mexErrMsgTxt("NS is not int32");}
	if (mxGetClassID(prhs[6]) != mxSINGLE_CLASS) {mexErrMsgTxt("rResol is not single");}

	// initialize host variables
	int const *NS = (int*)mxGetData(prhs[5]);			// total number of spheres
	float const *rResol = (float*)mxGetData(prhs[6]);	// resolution of lookups
	int const blocksPerGrid = (NS[0] + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
	
	// initialize device variables
	float const *d_real_x = (float const *)(mxGPUGetDataReadOnly(mx_real_x));
	float const *d_imag_x = (float const *)(mxGPUGetDataReadOnly(mx_imag_x));
	float const *d_real_h = (float const *)(mxGPUGetDataReadOnly(mx_real_h));
	float const *d_imag_h = (float const *)(mxGPUGetDataReadOnly(mx_imag_h));
	float const *d_sPos = (float const *)(mxGPUGetDataReadOnly(mx_sPos));
	float	 	*d_real_Wx = (float *)(mxGPUGetData(mx_real_Wx));
	float	 	*d_imag_Wx = (float *)(mxGPUGetData(mx_imag_Wx));
	
	// start computation
	for (int s2=1; s2<=NS[0]; s2++)
	{
		translationMatrixProduct<<< blocksPerGrid, THREADS_PER_BLOCK >>> (s2, NS[0], d_sPos,
																		  d_real_h, d_imag_h, rResol[0],
																		  d_real_x, d_imag_x, 
																		  d_real_Wx, d_imag_Wx);
																		
		//mexPrintf("sphere %i\n",s2);
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
	mxGPUDestroyGPUArray(mx_sPos);
	mxGPUDestroyGPUArray(mx_real_Wx);
	mxGPUDestroyGPUArray(mx_imag_Wx);
	
	//cudaDeviceReset();  // necessary for profiling
}
