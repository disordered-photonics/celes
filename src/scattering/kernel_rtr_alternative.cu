#include "wigner_D_CUDA.cuh"
#include "translation_coefficients_CUDA.cuh"


__global__ void translationMatrixProduct(int const s2, int const NS, float const *sPosArray,
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
