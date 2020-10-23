% plot MSTM output
clear; close all

% reload data from mat file if available, otherwise save it afterwards
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

    % read cut spheres positions
    spheres = dlmread(filename,'',[2, 0, cutspheres+1, 2]);
    
    % read near-field data for 0 deg polarization
    mstmOutput = dlmread(filename,'',cutspheres+2,0);
    xvec=unique(mstmOutput(:,1));
    zvec=unique(mstmOutput(:,2));
    E0 = cat(3, reshape(mstmOutput(:,3) + 1i*mstmOutput(:,4), [length(zvec), length(xvec)]), ...
                reshape(mstmOutput(:,5) + 1i*mstmOutput(:,6), [length(zvec), length(xvec)]), ...
                reshape(mstmOutput(:,7) + 1i*mstmOutput(:,8), [length(zvec), length(xvec)]));
    H0 = cat(3, reshape(mstmOutput(:,9) + 1i*mstmOutput(:,10),[length(zvec), length(xvec)]), ...
                reshape(mstmOutput(:,11)+ 1i*mstmOutput(:,12),[length(zvec), length(xvec)]), ...
                reshape(mstmOutput(:,13)+ 1i*mstmOutput(:,14),[length(zvec), length(xvec)]));
    E0 = single(E0); % the ascii output has only a few digits, therefore we can store it as single
    H0 = single(H0);
    
    % read near-field data for 90 deg polarization
    filename = 'mstm_benchmark_90deg_nf.dat';
    mstmOutput = dlmread(filename,'',cutspheres+2,0);
    E90= cat(3, rot90(reshape(mstmOutput(:,3) + 1i*mstmOutput(:,4), [length(zvec), length(xvec)]),0), ...
                rot90(reshape(mstmOutput(:,5) + 1i*mstmOutput(:,6), [length(zvec), length(xvec)]),0), ...
                rot90(reshape(mstmOutput(:,7) + 1i*mstmOutput(:,8), [length(zvec), length(xvec)]),0));
    H90= cat(3, rot90(reshape(mstmOutput(:,9) + 1i*mstmOutput(:,10),[length(zvec), length(xvec)]),0), ...
                rot90(reshape(mstmOutput(:,11)+ 1i*mstmOutput(:,12),[length(zvec), length(xvec)]),0), ...
                rot90(reshape(mstmOutput(:,13)+ 1i*mstmOutput(:,14),[length(zvec), length(xvec)]),0));
    E90 = single(E90);
    H90 = single(H90);
    [xgrid,zgrid]=meshgrid(xvec,zvec);
end

%%

% stack fields for convenience
fields0deg = cat(3, E0, sqrt(sum(abs(E0).^2, 3)), H0, sqrt(sum(abs(H0).^2, 3)));
fields90deg = cat(3, E90, sqrt(sum(abs(E90).^2, 3)), H90, sqrt(sum(abs(H90).^2, 3)));

% panel titles
comp = {'Re(E_x)', 'Re(E_y)', 'Re(E_z)', '|E|', ...
        'Re(H_x)', 'Re(H_y)', 'Re(H_z)', '|H|'};

% plot all components for 0deg polarization angle
figure('Name','MSTM, polarization angle = 0 deg', 'NumberTitle','off');
for sp=1:numel(comp)
    subplot(2,4,sp)
    custom_MSTM_panel(gca, xvec, zvec, real(fields0deg(:,:,sp)), spheres, comp{sp})
end
linkaxes(findall(gcf,'type','axes'))

% plot all components for 90deg polarization angle
figure('Name','MSTM, polarization angle = 90 deg', 'NumberTitle','off');
for sp=1:numel(comp)
    subplot(2,4,sp)
    custom_MSTM_panel(gca, xvec, zvec, real(fields90deg(:,:,sp)), spheres, comp{sp})
end
linkaxes(findall(gcf,'type','axes'))

if ~(exist('mstm_benchmark_nf.mat', 'file') == 2)
    clearvars cutspheres filename fld_mstm mstmOutput fields0deg fields90deg
    save('mstm_benchmark_nf.mat')
end

%%

function custom_MSTM_panel(ax, x, y, fld, spherecoords, titlestr)
% compose panel layout
imagesc(x,y,fld), axis equal tight
if contains(titlestr,'Re')
    cmap = interp1(1:3,[0 0 1; 1 1 1; 1 0 0],linspace(1,3,256));
    colormap(ax,cmap)
    caxis([-2, 2])
else
    caxis([0, 2])
end

pos = spherecoords(:,1:2);
r = spherecoords(:,3);
for s=1:size(r)
    rectangle(ax, ...
        'Position', [pos(s,:)-r(s), [2,2]*r(s)], ...
        'Curvature', [1 1], ...
        'FaceColor', 'none', ...
        'EdgeColor', [0,0,0], ...
        'LineWidth', 0.75)
end
title(titlestr)
xlabel('k*z'), ylabel('k*x')
view([270,90])
drawnow
end
