function fmm = assemble_neighborhood_coupling_matrices(fmm,ab5_table,positionArray,particleLmax,k)

for jbox1 = 1:fmm.numberOfBoxes  % box with receiving particles
    jbox1
    box1Positions = positionArray(fmm.particleIdcs{jbox1},:);
    NS1 = length(fmm.particleIdcs{jbox1});    
    for jbox2 = fmm.neighborBoxIdcs{jbox1} % box with emitting particles
        box2Positions = positionArray(fmm.particleIdcs{jbox2},:);
        NS2 = length(fmm.particleIdcs{jbox2});
        
        fmm.neighborhoodCouplingMatrices{jbox1}{jbox2} = zeros(jmult_max(length(fmm.particleIdcs{jbox1}),particleLmax),jmult_max(length(fmm.particleIdcs{jbox2}),particleLmax),'single');
        
        [x2,x1]=meshgrid(box2Positions(:,1),box1Positions(:,1));
        [y2,y1]=meshgrid(box2Positions(:,2),box1Positions(:,2));
        [z2,z1]=meshgrid(box2Positions(:,3),box1Positions(:,3));
        x1mnx2 = x1-x2;
        y1mny2 = y1-y2;
        z1mnz2 = z1-z2;
        dTab = sqrt(x1mnx2.^2+y1mny2.^2+z1mnz2.^2);
        ctTab = z1mnz2./dTab;
        stTab = sqrt(1-ctTab.^2);
        phiTab = atan2(y1mny2,x1mnx2);
        Plm = legendre_normalized_trigon(ctTab,stTab,2*particleLmax);
        
        for p=0:2*particleLmax
            sphHank = sph_bessel(3,p,k*dTab);
            for tau1=1:2
                for l1=1:particleLmax
                    for m1=-l1:l1
                        n1=multi2single_index(1,tau1,l1,m1,particleLmax);
                        n1S1Arr=(1:NS1)+(n1-1)*NS1;
                        for tau2=1:2
                            for l2=1:particleLmax
                                for m2=-l2:l2
                                    if abs(m1-m2)<=p
                                        n2=multi2single_index(1,tau2,l2,m2,particleLmax);
                                        n2S2Arr=(1:NS2)+(n2-1)*NS2;
                                        W = ab5_table(n2,n1,p+1) * Plm{p+1,abs(m1-m2)+1} .* sphHank .* exp(1i*(m2-m1)*phiTab) ;
                                        if jbox1==jbox2
                                            s1eqs2=logical(eye(NS1));
                                            W(s1eqs2(:))=0; % jS1=jS2
                                        end
                                        fmm.neighborhoodCouplingMatrices{jbox1}{jbox2}(n1S1Arr,n2S2Arr) = W;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end