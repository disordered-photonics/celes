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

%===============================================================================
%> @brief Show an image-plot of the near field
%>
%> The sphere positions are displayed as circles.
%>
%> @param   ax (axes object): axes to plot in
%>
%> @param   simulation (celes_simulation): simulation object
%>          containing the solution to plot
%>
%> @param   component (string): select 'real Ex', 'real Ey', 'real Ez'
%>          or 'abs E'
%>
%> @param   fieldType (string): select 'Total field', 'Scattered field'
%>          or 'Initial field'
%===============================================================================
function plot_field(ax,simulation,component,fieldType)

hold(ax,'on')

pArr = simulation.input.particles.positionArray;
rArr = simulation.input.particles.radiusArray;
dims = simulation.output.fieldPointsArrayDims;

switch lower(fieldType)
    case 'total field'
        E = simulation.output.totalField;
    case 'scattered field'
        E = simulation.output.scatteredField + simulation.output.internalField;
    case 'initial field'
        E = simulation.output.initialField;
end

switch lower(component)
    case 'real ex'
        fld = reshape(gather(real(E(:,1))), dims);
    case 'real ey'
        fld = reshape(gather(real(E(:,2))), dims);
    case 'real ez'
        fld = reshape(gather(real(E(:,3))), dims);
    case 'abs e'
        fld = reshape(gather(sqrt(sum(abs(E).^2,2))), dims);
end


fldPoints = reshape([simulation.output.fieldPoints(:,1), ...
                     simulation.output.fieldPoints(:,2), ...
                     simulation.output.fieldPoints(:,3)], [dims,3]);

if all(fldPoints(:,:,1) == fldPoints(1,1,1))    % fldPoints are on the yz plane
    perpdim = 1;                                % 1->x is the perp. direction
    draw_image(ax, fldPoints, fld, perpdim, pArr, rArr)
    xlabel('y')
    ylabel('z')
elseif all(fldPoints(:,:,2) == fldPoints(1,1,2))% fldPoints are on the xz plane
    perpdim = 2;                                % 2->y is the perp. direction
    draw_image(ax, fldPoints, fld, perpdim, pArr, rArr)
    xlabel('x')
    ylabel('z')
elseif all(fldPoints(:,:,3) == fldPoints(1,1,3))% fldPoints are on the xy plane
    perpdim = 3;                                % 3->z is the perp. direction
    draw_image(ax, fldPoints, fld, perpdim, pArr, rArr)
    xlabel('x')
    ylabel('y')
else
    error('fieldPoint must define an xy, xz or yz-like plane')
end

ax.DataAspectRatio = [1,1,1];
title([fieldType,', ',component])
hold(ax,'off')
end

function draw_image(ax, fldP, fld, perpdim, pArr, rArr)
xy = setdiff([1,2,3], perpdim);                 % here xy are the in-plane dimensions
x = fldP(:,:,xy(1));
y = fldP(:,:,xy(2));
imagesc(x(1,:), y(:,1), fld)                    % plot field on a xy plane
dist = abs(pArr(:,perpdim) - fldP(1,1,perpdim));% particle distances from xy plane
idx = find(dist<rArr);                          % find particles intersecting the plane
rArr(idx) = sqrt(rArr(idx).^2 - dist(idx).^2);  % overwrite radius of the intersection circle
for i=1:length(idx)
    rectangle(ax, ...
             'Position', [pArr(idx(i),xy)-rArr(idx(i)), [2,2]*rArr(idx(i))], ...
             'Curvature', [1 1], ...
             'FaceColor', 'none', ...
             'EdgeColor', [0,0,0], ...
             'LineWidth', 0.75)
end
axis([min(x(:)),max(x(:)),min(y(:)),max(y(:))]) % set axis tight to fldPoints
end