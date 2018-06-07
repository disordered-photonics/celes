/*=============================================================================
@brief 	Evaluate the lookup of the spherical Bessel or Neumann function with 
		cubic spline interpolation.
		
		The first two entries of the distances array should both be zero,
		like this:
		
		>> lookupParticleDistances = [0, 0:step:maxDistance];
		
		The lookup table is assumed to be of the following shape:
		
		>> besselTable = gpuArray.zeros(2*lmax+2, length(lookupParticleDistances), 'single');
		>> neumannTable = gpuArray.zeros(2*lmax+2, length(lookupParticleDistances), 'single');
		>> for p = 0:(2*lmax+1)
		>>	 spbs = sph_bessel(3,p,lookupParticleDistances);
		>>	 besselTable(p+1,:) = real(spbs);
		>>	 neumannTable(p+1,:) = imag(spbs);
		>> end

@param	p 			Spherical Bessel/Neumann function order
@param	r			Radial position
@param	spjTable	Pointer to lookup table
@param	rResol		Sampling distance of radial position

@retval 			Interpolated value
=============================================================================*/
__device__ float interpolateLookup(int const p, float const r, 
                                   float const *spjTable, float const rResol)
{
	float spj = 0.0f;
	float rPos = r/rResol;
	int rIdx = int(rPos);    					// points to table position -1, because for each p, the first entry with respect to r in the spjTable is copied 
	rPos -= rIdx; 							 	// (remainder of r/rResol) / rResol
	float rPos2 = pow(rPos,2);
	float rPos3 = pow(rPos,3);
	spj = ((-rPos3+2*rPos2-rPos) * spjTable[rIdx*(2*LMAX+2)+p]
			+ (3*rPos3-5*rPos2+2) * spjTable[(rIdx+1)*(2*LMAX+2)+p]
			+ (-3*rPos3+4*rPos2+rPos) * spjTable[(rIdx+2)*(2*LMAX+2)+p]
			+ (rPos3-rPos2) * spjTable[(rIdx+3)*(2*LMAX+2)+p])/2;
	return spj;
}


__device__ int idx(const int m, const int l, const int lprime)
{
	return m + l*(LMAX+1) + lprime*(LMAX+1)*(LMAX+1);
}

__device__ void axialTranslationCoefficients(const float kz,
                                             float const *jTable, float const *yTable,
											 const float rResol,
                                             float *real_A, float *imag_A, 
                                             float *real_B, float *imag_B)
{
	float real_C[(LMAX+1)*(2*LMAX+2)*(LMAX+1)];		// TODO: reduce register usage by more intelligent indexing
	float imag_C[(LMAX+1)*(2*LMAX+2)*(LMAX+1)];
	
	for (int l=0; l<=LMAX; l++){
		for (int lprime=0; lprime<=(2*LMAX+1); lprime++){
			for (int m=min(l, lprime); m<=LMAX; m++){	// TODO: check if this breaks constant indexing (important for register usage)
				real_C[idx(m, l, lprime)] = 0.0f;
				imag_C[idx(m, l, lprime)] = 0.0f;
			}
		}
	}
	
	// [Doicu Appendix B (page 280)]
	for (int lprime=0; lprime<=(2*LMAX+1); lprime++){
		real_C[idx(0, 0, lprime)] = powf(-1.0f, lprime) * sqrtf(2.0f*lprime+1.0f) * interpolateLookup(lprime, kz, jTable, rResol);
		imag_C[idx(0, 0, lprime)] = powf(-1.0f, lprime) * sqrtf(2.0f*lprime+1.0f) * interpolateLookup(lprime, kz, yTable, rResol);
	}
	
	for (int m=1; m<=LMAX; m++){
		for (int lprime=m; lprime<=(2*LMAX); lprime++){
			float const fac1 = sqrtf( (2.0f*m+1.0f) / (2.0f*m) );
			float const fac2 = sqrtf( (lprime+m-1.0f) * (lprime+m) / (2.0f*lprime-1.0f) / (2.0f*lprime+1.0f) );
			float const fac3 = sqrtf( (lprime-m+1.0f) * (lprime-m+2.0f) / (2.0f*lprime+1.0f) / (2.0f*lprime+3.0f) );
			real_C[idx(m, m, lprime)] = fac1 * (fac2 * real_C[idx(m-1, m-1, lprime-1)] + fac3 * real_C[idx(m-1, m-1, lprime+1)]);
			imag_C[idx(m, m, lprime)] = fac1 * (fac2 * imag_C[idx(m-1, m-1, lprime-1)] + fac3 * imag_C[idx(m-1, m-1, lprime+1)]);
		}
	}
    
	// [Doicu B.65 (page 279)]
	for (int m=0; m<=LMAX; m++){
		for (int l=m; l<=(LMAX-1); l++){  // we write the coefficient C_{m l+1, m lprime}
			for (int lprime=m; lprime<=(2*LMAX-l); lprime++){
				float real_term1;
				float imag_term1;
				if (l==m) {
					real_term1 = 0.0f; imag_term1 = 0.0f; 
				}
				else {
					real_term1 = sqrtf( (l-m)*(l+m) / ((2.0f*l-1.0f) * (2.0f*l+1.0f)) ) * real_C[idx(m, l-1, lprime)];
					imag_term1 = sqrtf( (l-m)*(l+m) / ((2.0f*l-1.0f) * (2.0f*l+1.0f)) ) * imag_C[idx(m, l-1, lprime)];
				}
				
				float real_term2;
				float imag_term2;
				if (lprime==m) {
					real_term2 = 0.0f; imag_term2 = 0.0f;
				}
				else {
					real_term2 = sqrtf( (lprime-m)*(lprime+m) / ((2.0f*lprime-1.0f) * (2.0f*lprime+1.0f)) ) * real_C[idx(m, l, lprime-1)];
					imag_term2 = sqrtf( (lprime-m)*(lprime+m) / ((2.0f*lprime-1.0f) * (2.0f*lprime+1.0f)) ) * imag_C[idx(m, l, lprime-1)];
				}
				
				float const real_term3 = - sqrtf( (lprime-m+1.0f)*(lprime+m+1.0f) / ((2.0f*lprime+1.0f) * (2.0f*lprime+3.0f)) ) 
				                           * real_C[idx(m, l, lprime+1)];
				float const imag_term3 = - sqrtf( (lprime-m+1.0f)*(lprime+m+1.0f) / ((2.0f*lprime+1.0f) * (2.0f*lprime+3.0f)) ) 
				                           * imag_C[idx(m, l, lprime+1)];										   
				
				real_C[idx(m, l+1, lprime)] = (real_term1 + real_term2 + real_term3) / sqrtf( (l-m+1.0f)*(l+m+1.0f) / ((2.0f*l+1.0f)*(2.0f*l+3.0f)) );
				imag_C[idx(m, l+1, lprime)] = (imag_term1 + imag_term2 + imag_term3) / sqrtf( (l-m+1.0f)*(l+m+1.0f) / ((2.0f*l+1.0f)*(2.0f*l+3.0f)) );
			}
		}
	}
	
	// [Doicu Appendix B (page 285)]
	for (int m=0; m<=LMAX; m++){
		for (int lprime=max(1, m); lprime<=LMAX; lprime++){
			for (int l=max(1, m); l<=LMAX; l++){
				float const prefac = sqrtf( lprime*(lprime+1.0f) / (l*(l+1.0f)) );
				
				float const real_term2 = kz/(lprime+1.0f) * sqrtf((lprime-m+1.0f)*(lprime+m+1.0f) / ((2.0f*lprime+1.0f)*(2.0f*lprime+3.0f))) * real_C[idx(m, l, lprime+1)];
				float const imag_term2 = kz/(lprime+1.0f) * sqrtf((lprime-m+1.0f)*(lprime+m+1.0f) / ((2.0f*lprime+1.0f)*(2.0f*lprime+3.0f))) * imag_C[idx(m, l, lprime+1)];
				
				float const real_term3 = kz/lprime * sqrtf((lprime-m)*(lprime+m) / ((2.0f*lprime+1.0f)*(2.0f*lprime-1.0f))) * real_C[idx(m, l, lprime-1)];
				float const imag_term3 = kz/lprime * sqrtf((lprime-m)*(lprime+m) / ((2.0f*lprime+1.0f)*(2.0f*lprime-1.0f))) * imag_C[idx(m, l, lprime-1)];
				
				real_A[idx(m,l,lprime)] = prefac * (real_C[idx(m,l,lprime)] + real_term2 + real_term3);
				imag_A[idx(m,l,lprime)] = prefac * (imag_C[idx(m,l,lprime)] + imag_term2 + imag_term3);
				
				real_B[idx(m,l,lprime)] = -kz * m / sqrtf(l*lprime*(l+1.0f)*(lprime+1.0f)) * imag_C[idx(m,l,lprime)];
				imag_B[idx(m,l,lprime)] = kz * m / sqrtf(l*lprime*(l+1.0f)*(lprime+1.0f)) * real_C[idx(m,l,lprime)];
			}
		}
	}
}
