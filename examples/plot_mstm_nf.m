% plot MSTM output
clear; close all

if exist('mstm_benchmark_nf.mat', 'file') == 2
    load('mstm_benchmark_nf.mat');
else
    % get the number of lines to skip (# spheres cut by the near-field plane)
    filename = 'mstm_benchmark_0deg_nf.dat';
    fileID = fopen(filename,'r');
    formatSpec = '%f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 1, 'HeaderLines', 1);
    cutspheres = dataArray{:, 1};
    fclose(fileID);
    clearvars formatSpec fileID dataArray ans;
    
    % read near-field data for 0 deg polarization
    mstmOutput = dlmread(filename,'',cutspheres+2,0);
    xvec=unique(mstmOutput(:,1));
    zvec=unique(mstmOutput(:,2));
    E0 = cat(3, rot90(reshape(mstmOutput(:,3) + 1i*mstmOutput(:,4), [length(zvec), length(xvec)])), ...
                rot90(reshape(mstmOutput(:,5) + 1i*mstmOutput(:,6), [length(zvec), length(xvec)])), ...
                rot90(reshape(mstmOutput(:,7) + 1i*mstmOutput(:,8), [length(zvec), length(xvec)])));
    H0 = cat(3, rot90(reshape(mstmOutput(:,9) + 1i*mstmOutput(:,10),[length(zvec), length(xvec)])), ...
                rot90(reshape(mstmOutput(:,11)+ 1i*mstmOutput(:,12),[length(zvec), length(xvec)])), ...
                rot90(reshape(mstmOutput(:,13)+ 1i*mstmOutput(:,14),[length(zvec), length(xvec)])));
    E0 = single(E0); % the ascii output has only a few digits...
    H0 = single(H0);
    
    % read near-field data for 90 deg polarization
    filename = 'mstm_benchmark_90deg_nf.dat';
    mstmOutput = dlmread(filename,'',cutspheres+2,0);
    E90= cat(3, rot90(reshape(mstmOutput(:,3) + 1i*mstmOutput(:,4), [length(zvec), length(xvec)])), ...
                rot90(reshape(mstmOutput(:,5) + 1i*mstmOutput(:,6), [length(zvec), length(xvec)])), ...
                rot90(reshape(mstmOutput(:,7) + 1i*mstmOutput(:,8), [length(zvec), length(xvec)])));
    H90= cat(3, rot90(reshape(mstmOutput(:,9) + 1i*mstmOutput(:,10),[length(zvec), length(xvec)])), ...
                rot90(reshape(mstmOutput(:,11)+ 1i*mstmOutput(:,12),[length(zvec), length(xvec)])), ...
                rot90(reshape(mstmOutput(:,13)+ 1i*mstmOutput(:,14),[length(zvec), length(xvec)])));
    E90 = single(E90);
    H90 = single(H90);
    [xgrid,zgrid]=meshgrid(xvec,zvec);
end

%%

cmap = interp1(1:3,[0 0 1; 1 1 1; 1 0 0],linspace(1,3,256)); % define a diverging colormap

figure
subplot(2,4,4)
imagesc(xvec, zvec, sqrt(sum(abs(E0).^2, 3))), axis equal tight
caxis([0,2])
title('|E|, 0 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,1)
imagesc(xvec, zvec, real(E0(:,:,1))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('E_x, 0 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,2)
imagesc(xvec, zvec, real(E0(:,:,2))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('E_y, 0 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,3)
imagesc(xvec, zvec, real(E0(:,:,3))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('E_z, 0 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,8)
imagesc(xvec, zvec, sqrt(sum(abs(H0).^2, 3))), axis equal tight
caxis([0,2])
title('|H|, 0 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,5)
imagesc(xvec, zvec, real(H0(:,:,1))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('H_x, 0 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,6)
imagesc(xvec, zvec, real(H0(:,:,2))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('H_y, 0 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,7)
imagesc(xvec, zvec, real(H0(:,:,3))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('H_z, 0 deg'), xlabel('k*x'), ylabel('k*z')
linkaxes(findall(gcf,'type','axes'))

figure
subplot(2,4,4)
imagesc(xvec, zvec, sqrt(sum(abs(E90).^2, 3))), axis equal tight
caxis([0,2])
title('|E|, 90 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,1)
imagesc(xvec, zvec, real(E90(:,:,1))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('E_x, 90 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,2)
imagesc(xvec, zvec, real(E90(:,:,2))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('E_y, 90 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,3)
imagesc(xvec, zvec, real(E90(:,:,3))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('E_z, 90 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,8)
imagesc(xvec, zvec, sqrt(sum(abs(H90).^2, 3))), axis equal tight
caxis([0,2])
title('|H|, 90 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,5)
imagesc(xvec, zvec, real(H90(:,:,1))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('H_x, 90 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,6)
imagesc(xvec, zvec, real(H90(:,:,2))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('H_y, 90 deg'), xlabel('k*x'), ylabel('k*z')

subplot(2,4,7)
imagesc(xvec, zvec, real(H90(:,:,3))), axis equal tight
caxis([-2,2])
colormap(gca,cmap)
title('H_z, 90 deg'), xlabel('k*x'), ylabel('k*z')
linkaxes(findall(gcf,'type','axes'))

clearvars cutspheres filename fld_mstm mstmOutput;
if ~(exist('mstm_benchmark_nf.mat', 'file') == 2)
    save('mstm_benchmark_nf.mat')
end