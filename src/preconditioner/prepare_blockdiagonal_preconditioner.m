%  Copyright (c) 2017, Amos Egel (KIT), Lorenzo Pattelli (LENS)
%                      Giacomo Mazzamuto (LENS)
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%
%  * Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%
%  * Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
%  * Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

% %======================================================================
% %> @brief Compute the master blockmatrix for each partition and store it as
% %> well as its LU decomposition
% %>
% %> @param       simulation(celes_simulation)
% %>
% %> @retval      simulation(celes_simulation): simulation object with
% %>              updated
% %>              fall into the same partition cuboid
% %======================================================================
% function simulation = prepare_blockdiagonal_preconditioner(simulation)
% 
% fprintf(1,'prepare blockdiagonal preconditioner ...\n');
% msg = '';
% lmax=simulation.numerics.lmax;
% nmax=simulation.numerics.nmax;
% k = simulation.input.k_medium;
% partitioning=simulation.numerics.solver.preconditioner.partitioning;
% 
% for jp=1:length(partitioning)
%     spherArr=partitioning{jp};
%     NSi = length(spherArr);
%     
%     Idcs= [];
%     for n=1:nmax
%         Idcs = [Idcs;spherArr+simulation.input.particles.number*(n-1)];
%     end
%     simulation.numerics.solver.preconditioner.partitioningIdcs{jp} = Idcs;
%     
%     fprintf(1,repmat('\b',[1,length(msg)]));
%     msg = sprintf('partition %i of %i: %i particles -- compute master matrix ...',jp,length(partitioning),length(spherArr));
%     fprintf(1,msg);
%     
%     M = eye(NSi*simulation.numerics.nmax,'single');
%     
%     [x2,x1]=meshgrid(simulation.input.particles.positionArray(spherArr,1));
%     [y2,y1]=meshgrid(simulation.input.particles.positionArray(spherArr,2));
%     [z2,z1]=meshgrid(simulation.input.particles.positionArray(spherArr,3));
%     x1mnx2 = x1-x2;
%     y1mny2 = y1-y2;
%     z1mnz2 = z1-z2;
%     dTab = sqrt(x1mnx2.^2+y1mny2.^2+z1mnz2.^2);    
%     ctTab = z1mnz2./dTab;
%     stTab = sqrt(1-ctTab.^2);
%     phiTab = atan2(y1mny2,x1mnx2);
%     Plm = legendre_normalized_trigon(ctTab,stTab,2*lmax);
%     
%     for p=0:2*lmax
%         sphHank = sph_bessel(3,p,k*dTab);
%         for tau1=1:2
%             for l1=1:lmax
%                 for m1=-l1:l1
%                     n1=multi2single_index(1,tau1,l1,m1,lmax);
%                     n1S1Arr=(1:NSi)+(n1-1)*NSi;
%                     for tau2=1:2
%                         for l2=1:lmax
%                             for m2=-l2:l2
%                                 if abs(m1-m2)<=p
%                                     n2=multi2single_index(1,tau2,l2,m2,lmax);
%                                     n2S2Arr=(1:NSi)+(n2-1)*NSi;
%                                     TWn1n2 = simulation.tables.mieCoefficients(n1)*simulation.tables.translationTable.ab5(n2,n1,p+1) * Plm{p+1,abs(m1-m2)+1} .* sphHank .* exp(1i*(m2-m1)*phiTab) ;
%                                     s1eqs2=logical(eye(NSi));
%                                     TWn1n2(s1eqs2(:))=0; % jS1=jS2
%                                     M(n1S1Arr,n2S2Arr) = M(n1S1Arr,n2S2Arr) - TWn1n2;
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     
%     fprintf(1,repmat('\b',[1,length(msg)]));
%     msg = sprintf('partition %i of %i: %i particles -- factorize master matrix ...',jp,length(partitioning),length(spherArr));
%     fprintf(1,msg);
%     
%     [Y,U,P] = lu(simulation.numerics.deviceArray(M),'vector');
%     Y(U~=0) = U(U~=0); 
%     simulation.numerics.solver.preconditioner.masterMatrices{jp}=M;
%     simulation.numerics.solver.preconditioner.factorizedMasterMatrices{jp}.Y=gather(Y);
%     simulation.numerics.solver.preconditioner.factorizedMasterMatrices{jp}.P=gather(P);
% end
% fprintf(' done\n')