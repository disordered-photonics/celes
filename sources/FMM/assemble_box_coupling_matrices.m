function fmm = assemble_box_coupling_matrices(fmm,ab5_table,k)


relativePositions = fmm.boxRelativeOffsetIdcs * fmm.edgeSize;

fmm.boxCouplingMatrices = zeros(jmult_max(1,fmm.lmax),jmult_max(1,fmm.lmax),length(relativePositions(:,1)),'single');

d = sqrt(sum(relativePositions.^2,2));
ct = relativePositions(:,3)./d;
st = sqrt(1-ct.^2);
phi = atan2(relativePositions(:,2),relativePositions(:,1));
Plm = legendre_normalized_trigon(ct,st,2*fmm.lmax);

for p=0:(2*fmm.lmax)
    sphBess = sph_bessel(3,p,k*d);
    for tau1=1:2
        for l1=1:fmm.lmax
            for m1=-l1:l1
                n1=multi2single_index(1,tau1,l1,m1,fmm.lmax);
                for tau2=1:2
                    for l2=1:fmm.lmax
                        for m2=-l2:l2
                            if abs(m1-m2)<=p
                                n2=multi2single_index(1,tau2,l2,m2,fmm.lmax);
                                fmm.boxCouplingMatrices(n1,n2,:) = ab5_table(n2,n1,p+1) * Plm{p+1,abs(m1-m2)+1} .* sphBess .* exp(1i*(m2-m1)*phi) ;
                            end
                        end
                    end
                end
            end
        end
    end
end

