% function T = T_matrix(model)
% 
% 
% T=zeros(1,jmult_max(1,model.lmax),'single');
% 
% for tau=1:2
%     for l=1:model.lmax
%         for m=-l:l
%             jmult = multi2single_index(1,tau,l,m,model.lmax);
%             T(jmult) = T_entry(tau,l,model.medium_k,model.particle_k,model.particle_radius);
%         end
%     end
% end