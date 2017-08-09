function T = T_block_NFDS_spheroid(wavelength,nC,nM,C,AB,lmax,Nrank,nInt)
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

%nInt=200;
useDS=true;
% useDS=false;
namestring='spheroid';

if isunix
    cd scattering/NFM-DS_linux/INPUTFILES/
    write_NFMDS_input_spheroid(wavelength,nM,nC,C,AB,nInt,useDS,namestring,Nrank)
    cd ../TMATSOURCES/
    dlmwrite('one.txt',1)
    system('./1_tmatrix_NFDS<one.txt')
    %system('./1_tmatrix_NFDS')
    cd ../TMATFILES/
    TCalc = read_NFMDS_output('TmatForTSPLspheroid.dat',Nrank);
    cd ../../..
else
    cd scattering/NFM-DS/INPUTFILES/
    write_NFMDS_input_spheroid(wavelength,nM,nC,C,AB,nInt,useDS,namestring,Nrank)
    cd ../TMATSOURCES/
    dlmwrite('one.txt',1)
    system('1_T_Matrix_Fortran_Program.exe<one.txt')
    cd ../TMATFILES/
    TCalc = read_NFMDS_output('TmatForTSPLspheroid.dat',Nrank);
    cd ../../..
end

nmax=jmult_max(1,lmax);
T=zeros(nmax);

for tau=1:2
    for l=1:lmax
        for m=-l:l
            n1 = multi2single_index(1,tau,l,m,lmax);
            n1Calc = multi2single_index(1,tau,l,m,Nrank);
            for tau2=1:2
                for l2=1:lmax
                    for m2=-l2:l2
                        n2 = multi2single_index(1,tau2,l2,m2,lmax);
                        n2Calc = multi2single_index(1,tau2,l2,m2,Nrank);
                        T(n1,n2)=TCalc(n1Calc,n2Calc);
                    end
                end
            end
        end
    end
end



