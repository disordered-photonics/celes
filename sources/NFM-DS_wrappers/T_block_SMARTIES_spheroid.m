function T = T_block_SMARTIES_spheroid(wavelength,nC,nM,C,AB,lmax,Nrank,nInt)
%
% wavelength: vacuum wavelength
% nC: complex refractive index of cylinder
% nM: refractive index of medium
% C: half axis in z-direction
% AB: half axis in xy-direction
% lmax: maximal multipole order
% Nrank: maximial multipole order for NFMDS algorithm 
% nInt: numer of integration points for NFMDS algorithm

% output
% T: matrix of dimension nind x nind, where nind is the total number of
% indexcombinations

stParams.a=AB;
stParams.c=C;
stParams.k1=2*pi/wavelength*nM;
stParams.s=nC/nM;


%% Convergence parameters
% N:        Maximum multipole order for T-matrix and series expansions of fields
% nNbTheta: Number of points for Gaussian quadratures to compute integrals in P and Q matrices

% Those can be estimated automatically for some desired accuracy as follows
% [N, nNbTheta] = sphEstimateNandNT(stParams, stOptions, 1e-8);

% In many instances, it will be more efficient to set those manually, e.g.
% N = 30;
% nNbTheta = 120;

% Add those to the parameters structure
stParams.N=lmax; 
stParams.nNbTheta=nInt;
stOptions = [];

[stCoa, CstTRa] = slvForT(stParams, stOptions);

for jj=1:(lmax+1)
    [Tm{jj},nvecT]=rvhGetFullMatrix(CstTRa{jj},'st4MT');
end

Tsmt = cell2mat(full_Tmat_SMARTIES(Tm,lmax));

%lookups.index = index_table_dominik(lmax);
T = reorganize_NFMDS_nonaxsym_TMat(Tsmt,lmax);
%T = TMat_SMARTIES_reorganization(Tsmt,lookups);