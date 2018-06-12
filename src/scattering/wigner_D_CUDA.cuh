/*=============================================================================
@brief 	Legendre polynomial by recursion

@param	x 			Argument of polynomial
@param	P			Float array of length LMAX+1 into which result is written
=============================================================================*/
__device__ void legendrePolynomial(float const x, float *P) 
{
	P[0] = 1;
	P[1] = x;
	for (int l=1; l<=(LMAX-1); l++)
	{
		P[l+1] = ((2*l + 1) * x * P[l] - l * P[l-1]) / (l+1);
	}
}

/*=============================================================================
@brief 	Quotient of two factorials n!/m!
		This function can be used to compute the ordinary factorial by 
		setting some m with 1 < m < 2

@param	n 			Argument of polynomial in the numerator
@param	m			Argument of polynomial in denominator

@retval n!/n2! as float
=============================================================================*/
__device__ float truncatedFactorial(float n, float const n2)
{
	float result = 1.0f;
	for (;n>=n2;n--){
		result *= n;
	}
	return result;
}

// NOTE: the following recursive factorial leads to memory errors with lmax>2
// ... maybe due stack size??
__device__ float truncatedFactorial_recursive(float const n, float const n2)
{
	if (n < n2) return 1.0f; else return n*truncatedFactorial(n-1.0f, n2);
}

/*=============================================================================
@brief 	SVWF rotation coefficients (Wigner-D symbols)
		The calculation goes along the following textbooks:
		[Doicu] 		Doicu, Adrian, Thomas Wriedt, and Yuri A. Eremin. 
						"Light scattering by systems of particles: null-field 
						method with discrete sources: theory and programs." 
						Springer, 2006.
	    [Mishchenko]	Mishchenko, Michael I., Larry D. Travis, 
						and Andrew A. Lacis. "Scattering, absorption, and 
						emission of light by  small particles." 
						Cambridge university press, 2002.

@param	m			azimuthal quantum number in rotated coordinate system
@param	m_prime		azimuthal quantum number in original coordinate system
@param	alpha		alpha Euler angle
@param	cosBeta		cosine of beta Euler angle
@param	gamma		gamma Euler angle
@param	D			Float array of length LMAX+1 into which result is written
					D[n] contains the symbol D^n_{m,m'}(\alpha,\beta,\gamma)
					for n=0...LMAX
=============================================================================*/
__device__ void wignerD(int const m, int const m_prime, 
						float const alpha, float const cosBeta, float const gamma, 
						float *real_D, float *imag_D)
{
	// write Wigner-d symbols into real_D
	if ((m==0) && (m_prime==0)) {
		legendrePolynomial(cosBeta, real_D);
	} 
	else {
		int const smin = max(abs(m), abs(m_prime));
		
		float zeta;
		if (m_prime >= m) zeta = 1.0f; else zeta = powf(-1.0f, m-m_prime);  // [Mishchenko, B.16 (page 364)]
		
		float maxm;
		float minm;
		if (abs(m-m_prime) > abs(m+m_prime)) {
			maxm = (float) abs(m-m_prime);
			minm = (float) abs(m+m_prime);
		} 
		else {
			maxm = (float) abs(m+m_prime);
			minm = (float) abs(m-m_prime);
		}
	
		for (int s=0; s<smin; s++) real_D[s] = 0.0f;  // this line can be removed
		// [Mishchenko, B.24 (page 365)]
		
		real_D[smin] = zeta * powf(2.0f, -smin)
		               * sqrtf( truncatedFactorial((float) 2*smin, maxm+0.5f) / truncatedFactorial(minm, 1.5f) )
					   //* sqrtf( (float) 2*smin / 1.5f)
					   * powf(1.0f - cosBeta, abs(m-m_prime)/2.0f) 
					   * powf(1.0f + cosBeta, abs(m+m_prime)/2.0f);
		
		
		// [Mishchenko, B.22 (page 365)]
		if (smin < LMAX){
			real_D[smin+1] = (2.0f*smin+1.0f) * (smin*(smin+1.0f) * cosBeta - m * m_prime) * real_D[smin]
                              / (smin * sqrtf((smin+1)*(smin+1)-m*m) * sqrtf((smin+1)*(smin+1) - m_prime*m_prime));	
		}
		
		for (int s=smin+1; s<=(LMAX-1); s++) {
			real_D[s+1] = ((2.0f*s+1.0f) * (s*(s+1.0f) * cosBeta - m * m_prime) * real_D[s]
			                - (s+1.0f) * sqrtf(s*s - m*m) * sqrtf(s*s - m_prime*m_prime) * real_D[s-1])
					       / (s * sqrtf((s+1)*(s+1)-m*m) * sqrtf((s+1)*(s+1) - m_prime*m_prime));
		}
		
	}
	
	// [Doicu, B.36 (page 271)]
	float delta = 1.0f;
	if ((m<0) && (m % 2)) delta *= -1.0f;
	if ((m_prime<0) && (m_prime % 2)) delta *= -1.0f;

	// [Doicu, B.34 (page 271)]
	for (int s=0; s<=LMAX; s++) {
		imag_D[s] = powf(-1.0f, m+m_prime) * delta * sinf(m*alpha + m_prime*gamma) * real_D[s];
		real_D[s] = powf(-1.0f, m+m_prime) * delta * cosf(m*alpha + m_prime*gamma) * real_D[s];
	}
}
