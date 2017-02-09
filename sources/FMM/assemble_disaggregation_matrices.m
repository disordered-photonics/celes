function fmm = assemble_disaggregation_matrices(fmm,positionArray,particleLmax,ab5_table,k)

for jbox = 1:fmm.numberOfBoxes
    NS = length(fmm.particleIdcs{jbox});
    fmm.particleDisaggregationMatrices{jbox} = zeros(jmult_max(1,particleLmax)*NS,jmult_max(1,fmm.lmax));
    relativePositionArray = bsxfun(@plus,positionArray(fmm.particleIdcs{jbox},:),-fmm.centerPositions{jbox});
    d = sqrt(sum(relativePositionArray.^2,2));
    ct = relativePositionArray(:,3)./d;
    st = sqrt(1-ct.^2);
    phi = atan2(relativePositionArray(:,2),relativePositionArray(:,1));
    Plm = legendre_normalized_trigon(ct,st,particleLmax+fmm.lmax);
    
    for p=0:(particleLmax+fmm.lmax)
        sphBess = sph_bessel(1,p,k*d);
        for tau1=1:2
            for l1=1:particleLmax
                for m1=-l1:l1
                    n1=multi2single_index(1,tau1,l1,m1,particleLmax);
                    n1SArr=(1:NS)+(n1-1)*NS;
                    for tau2=1:2
                        for l2=1:fmm.lmax
                            for m2=-l2:l2
                                if abs(m1-m2)<=p
                                    n2=multi2single_index(1,tau2,l2,m2,fmm.lmax);
                                    fmm.particleDisaggregationMatrices{jbox}(n1SArr,n2) = ab5_table(n2,n1,p+1) * Plm{p+1,abs(m1-m2)+1} .* sphBess .* exp(1i*(m2-m1)*phi) ;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
