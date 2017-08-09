function T = read_NFMDS_output(fileName,lmax)
% T = read_NFMDS_output(fileName,lmax)
% Reads the T-Matrix of the output of the NFM-DS software package,
% provided by Adrian Doicu and Thomas Wriedt  with the text book
% "Light Scattering by Systems of Particles"

infoData = fileread(['Info',fileName]);

if strfind(infoData,'The scatterer is an axisymmetric particle')  % is it output from TAXYM.f90?
    NrankPos = strfind(infoData,'Nrank =');
    Nrank = sscanf(infoData((NrankPos+7):end),'%i');
    MrankPos = strfind(infoData,'Mrank =');
    Mrank = sscanf(infoData((MrankPos+7):end),'%i');
    if Mrank<lmax
        warning('Mrank smaller than lmax')
    end
    TRaw = dlmread(fileName,'',3,0);
    sz = size(TRaw);
    for j1=1:(sz(2)/2)
        TNFMDS(:,j1) = TRaw(:,2*j1-1) + 1i*TRaw(:,2*j1);   % put two real numbers to one complex
    end
    if 2*Nrank > 10
        rowsize = ceil(2*Nrank/10)*10;
        numEl = length(TNFMDS(:));  % number of elements in raw data
        TNFMDS = reshape(TNFMDS.',rowsize,numEl/rowsize).';
        TNFMDS = TNFMDS(:,1:2*Nrank);
    end
    
    for m=-lmax:lmax
        Nmax_NFMDS = Nrank - max(1,abs(m)) + 1;
        for tau1=1:2
            for l1=max(1,abs(m)):lmax
                n1 = multi2single_index(1,tau1,l1,m,lmax);
                l1_NFMDS = l1 - max(1,abs(m)) + 1;
                n1_NFMDS = 2*Nrank*(abs(m)) + (tau1-1)*Nmax_NFMDS + l1_NFMDS;
                for tau2=1:2
                    for l2=max(1,abs(m)):lmax
                        n2 = multi2single_index(1,tau2,l2,m,lmax);
                        l2_NFMDS = l2 - max(1,abs(m)) + 1;
                        n2_NFMDS = (tau2-1)*Nmax_NFMDS + l2_NFMDS;
                        if abs(m)<=Mrank
                            if m>=0
                                T(n1,n2) = TNFMDS(n1_NFMDS,n2_NFMDS);
                            else
                                T(n1,n2) = TNFMDS(n1_NFMDS,n2_NFMDS) * (-1)^(tau1+tau2);
                            end
                        else
                            T(n1,n2) = 0;
                        end
                    end
                end
            end
        end
    end
else
    error('not yet implemented')
end
