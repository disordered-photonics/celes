#include "wigner_D_CUDA.cuh"
#include "translation_coefficients_CUDA.cuh"


__device__ void inner_loop(int const m, int const l2, 
                           int const s1, int const s2, int const NS,
						   float const phi, float const cosTheta,
						   float const kz,
						   float const *real_A, float const *imag_A, 
						   float const *real_B, float const *imag_B,
						   float const *real_x, float const *imag_x,
						   float *real_Wx, float *imag_Wx)
{	
	for (int m2=-l2; m2<=l2; m2++){
		float real_D1_array[LMAX+1];
		float imag_D1_array[LMAX+1];
		wignerD(m2, m, phi, cosTheta, 0, real_D1_array, imag_D1_array);
		float real_D1 = real_D1_array[l2];
		float imag_D1 = imag_D1_array[l2];									
	
		for (int m1=-LMAX; m1<=LMAX; m1++){
			float real_D2_array[LMAX+1];
			float imag_D2_array[LMAX+1];
			wignerD(m, m1, 0, cosTheta, -phi, real_D2_array, imag_D2_array);
			
			for (int l1=max(1, max(abs(m1),abs(m))); l1<=LMAX; l1++){
				float real_D2 = real_D2_array[l1] * powf(-1.0f, m+m1);  // [Doicu, B.41]
				float imag_D2 = imag_D2_array[l1] * powf(-1.0f, m+m1);
				
				for (int tau1=1; tau1<=2; tau1++){
					int n1 = (tau1-1)*LMAX*(LMAX+2)+(l1-1)*(l1+1)+l1+1+m1;
					int WxIdx = (n1-1)*NS+s1-1;
					float re_incr = 0.0f;
					float im_incr = 0.0f;
				
					for (int tau2=1; tau2<=2; tau2++){
						int n2 = (tau2-1)*LMAX*(LMAX+2)+(l2-1)*(l2+1)+l2+1+m2;
						int xIdx = (n2-1)*NS+s2-1;
						float re_xTmp = real_x[xIdx];
						float im_xTmp = imag_x[xIdx];
				
						float real_AB;
						float imag_AB;
						
						if (tau1==tau2){
							real_AB = real_A[l1];
							imag_AB = imag_A[l1];
						}
						else if (m<0){
							real_AB = -real_B[l1];
							imag_AB = -imag_B[l1];
						} 
						else {
							real_AB = real_B[l1];
							imag_AB = imag_B[l1];
						}
						
						float real_D1AB = real_D1*real_AB - imag_D1*imag_AB;
						float imag_D1AB = real_D1*imag_AB + imag_D1*real_AB;
						
						float real_D1ABD2 = real_D1AB * real_D2 - imag_D1AB * imag_D2;
						float imag_D1ABD2 = real_D1AB * imag_D2 + imag_D1AB * real_D2;
						
						re_incr += real_D1ABD2 * re_xTmp - imag_D1ABD2 * im_xTmp; 
						im_incr += real_D1ABD2 * im_xTmp + imag_D1ABD2 * re_xTmp;
					}
					
					real_Wx[WxIdx] += re_incr;
					imag_Wx[WxIdx] += im_incr;
				}
			}
		}
	}
}


__global__ void translationMatrixProduct(int const s2, int const NS, float const *sPosArray,
										 float const *sphericalBesselTable, float const *sphericalNeumannTable, float rResol,
										 float const *real_x, float const *imag_x, float *real_Wx, float *imag_Wx)
{
	//NOTE: dimensionless positions, i.e. k*xyz!!!
	int const s1 = blockDim.x * blockIdx.x + threadIdx.x + 1; // receiving sphere number (1...NS)
	float x21, y21, z21;
	float r, cosTheta, phi;
	
	if ((s1!=s2)&&s1<=NS)
	{
		// relative position
		x21 = sPosArray[3*(s1-1)]-sPosArray[3*(s2-1)];
		y21 = sPosArray[3*(s1-1)+1]-sPosArray[3*(s2-1)+1];
		z21 = sPosArray[3*(s1-1)+2]-sPosArray[3*(s2-1)+2];
		r = sqrtf(x21*x21+y21*y21+z21*z21);
		cosTheta = z21/r;
		phi = atan2f(y21,x21);

		// matrix vector product W*x -- remember that W is the transpose of A, such that indices are interchanged (l1<->l2 etc.)
				
		float real_A[LMAX+1] = {0.0f};  // SVWF axial translation coefficients, idx: l' (corresponding to l1)
		float imag_A[LMAX+1] = {0.0f};
		float real_B[LMAX+1] = {0.0f};
		float imag_B[LMAX+1] = {0.0f};
		float real_C[2*LMAX+2] = {0.0f};  // scalar SWF axial translation coefficients, idx: l' (corresponding to l1)
		float imag_C[2*LMAX+2] = {0.0f};  
		float real_C_mminus1[2*LMAX+2] = {0.0f};  // C_{m-1 l, m-1 l'}
		float imag_C_mminus1[2*LMAX+2] = {0.0f};  
		float real_C_lminus1[2*LMAX+2] = {0.0f};  // C_{m l-1, m l'}
		float imag_C_lminus1[2*LMAX+2] = {0.0f};
		float real_C_lminus2[2*LMAX+2] = {0.0f};  // C_{m l-2, m l'}
		float imag_C_lminus2[2*LMAX+2] = {0.0f};
		
		// treat m=0, l=0 in advance
		C00lprime(r, sphericalBesselTable, sphericalNeumannTable, rResol, real_C, imag_C);
		copy(real_C, real_C_mminus1);
		copy(imag_C, imag_C_mminus1);
		
		// the loop order is: m, l2, m2, m1, l1, tau1, tau2
		// reason for this twisted order:
		// m and l2 are needed for the iterative calculation of A,B
		// m and m2 are needed for the calculation of D1
		// m and m1 are needed for the calculation of D2
		for (int absm=0; absm<=LMAX; absm++){
			for (int l2=max(1,absm); l2<=LMAX; l2++){
				// compute axial translation coefficients A, B
				// l2 corresponds to l, whereas l1 corresponds to lprime
				if (l2==absm){
					Cmmlprime(absm, real_C_mminus1, imag_C_mminus1, real_C, imag_C);
					copy(real_C, real_C_mminus1);
					copy(imag_C, imag_C_mminus1);
					AB(absm, absm, r, real_C, imag_C, real_A, imag_A, real_B, imag_B);
					for (int lprime=0; lprime<2*LMAX+2; lprime++){ // reset C_{m l-1, m l'}
						real_C_lminus1[lprime] = 0.0f;
						imag_C_lminus1[lprime] = 0.0f;
					}
				} else {  // l2>absm
					copy(real_C_lminus1, real_C_lminus2);
					copy(imag_C_lminus1, imag_C_lminus2);
					copy(real_C, real_C_lminus1);
					copy(imag_C, imag_C_lminus1);
					Cmllprime(absm, l2, 
							  real_C_lminus2, imag_C_lminus2,
							  real_C_lminus1, imag_C_lminus1, 
							  real_C, imag_C);
					AB(absm, l2, r, real_C, imag_C, real_A, imag_A, real_B, imag_B);
				}
				// end of A, B computation
				
				if (absm==0) inner_loop(absm, l2, s1, s2, NS, phi, cosTheta, r, 
				                        real_A, imag_A, real_B, imag_B, 
										real_x, imag_x, real_Wx, imag_Wx);
				
				else {
					for (int m=-absm; m<=absm; m+=2*absm){
						inner_loop(m, l2, s1, s2, NS, phi, cosTheta, r, 
						           real_A, imag_A, real_B, imag_B, 
								   real_x, imag_x, real_Wx, imag_Wx);
					}
				}
			}
		}
	}
}


__global__ void translationMatrixProduct_precomputeAB(int const s2, int const NS, float const *sPosArray,
										              float const *sphericalBesselTable, float const *sphericalNeumannTable, float rResol,
													  float const *re_x, float const *im_x, float *re_Wx, float *im_Wx)
{
	//NOTE: dimensionless positions, i.e. k*xyz!!!
	int const s1 = blockDim.x * blockIdx.x + threadIdx.x + 1; // receiving sphere number (1...NS)
	float x21, y21, z21;
	float r, cosTheta, phi;
	
	float real_A[(LMAX+1)*(LMAX+1)*(LMAX+1)];  // TODO: more economic indexing
	float imag_A[(LMAX+1)*(LMAX+1)*(LMAX+1)];  
	float real_B[(LMAX+1)*(LMAX+1)*(LMAX+1)];  
	float imag_B[(LMAX+1)*(LMAX+1)*(LMAX+1)];  
	
	if ((s1!=s2)&&s1<=NS)
	{
		// relative position
		x21 = sPosArray[3*(s1-1)]-sPosArray[3*(s2-1)];
		y21 = sPosArray[3*(s1-1)+1]-sPosArray[3*(s2-1)+1];
		z21 = sPosArray[3*(s1-1)+2]-sPosArray[3*(s2-1)+2];
		
		r = sqrtf(x21*x21+y21*y21+z21*z21);
		cosTheta = z21/r;
		phi = atan2f(y21,x21);

		axialTranslationCoefficients(r, sphericalBesselTable, sphericalNeumannTable,
									 rResol, real_A, imag_A, real_B, imag_B);
		
		// matrix vector product W*x -- remember that W is the transpose of A (l1<->l2 etc.)
		for (int m=-LMAX; m<=LMAX; m++){
			for (int m1=-LMAX; m1<=LMAX; m1++){
				float real_D2_array[LMAX+1];
				float imag_D2_array[LMAX+1];
				wignerD(m, m1, 0, cosTheta, -phi, real_D2_array, imag_D2_array);
				
				for (int m2=-LMAX; m2<=LMAX; m2++){
					float real_D1_array[LMAX+1];
					float imag_D1_array[LMAX+1];
					wignerD(m2, m, phi, cosTheta, 0, real_D1_array, imag_D1_array);
					
					for (int l1=max(1, max(abs(m1),abs(m))); l1<=LMAX; l1++){
						float real_D2 = real_D2_array[l1] * powf(-1.0f, m+m1);  // [Doicu, B.41]
						float imag_D2 = imag_D2_array[l1] * powf(-1.0f, m+m1);

						for (int tau1=1; tau1<=2; tau1++){
							int n1 = (tau1-1)*LMAX*(LMAX+2)+(l1-1)*(l1+1)+l1+1+m1;
							int WxIdx = (n1-1)*NS+s1-1;
							float re_incr = 0.0f;
							float im_incr = 0.0f;
							
							for (int l2=max(1, max(abs(m2),abs(m))); l2<=LMAX; l2++){
								float real_D1 = real_D1_array[l2];
								float imag_D1 = imag_D1_array[l2];									
					
								for (int tau2=1; tau2<=2; tau2++){
									int n2 = (tau2-1)*LMAX*(LMAX+2)+(l2-1)*(l2+1)+l2+1+m2;
									int xIdx = (n2-1)*NS+s2-1;
									float re_xTmp = re_x[xIdx];
									float im_xTmp = im_x[xIdx];
							
									// axial translation
									float real_AB;
									float imag_AB;
									
									if (tau1==tau2){
										real_AB = real_A[idx(abs(m),l2,l1)];
										imag_AB = imag_A[idx(abs(m),l2,l1)];
									}
									else if (m<0){
										real_AB = -real_B[idx(-m,l2,l1)];
										imag_AB = -imag_B[idx(-m,l2,l1)];
									} 
									else {
										real_AB = real_B[idx(m, l2, l1)];
										imag_AB = imag_B[idx(m, l2, l1)];
									}
									
									float real_D1AB = real_D1*real_AB - imag_D1*imag_AB;
									float imag_D1AB = real_D1*imag_AB + imag_D1*real_AB;
									
									float real_D1ABD2 = real_D1AB * real_D2 - imag_D1AB * imag_D2;
									float imag_D1ABD2 = real_D1AB * imag_D2 + imag_D1AB * real_D2;
									
									re_incr += real_D1ABD2 * re_xTmp - imag_D1ABD2 * im_xTmp; 
									im_incr += real_D1ABD2 * im_xTmp + imag_D1ABD2 * re_xTmp;
								}
									//if (n2==1){
									//	re_Wx[WxIdx] = real_D2;
									//	im_Wx[WxIdx] = imag_D2;
									//}
							}
							re_Wx[WxIdx] += re_incr;
							im_Wx[WxIdx] += im_incr;
						}
					}
				} 
			}  
		}
	} 
}
