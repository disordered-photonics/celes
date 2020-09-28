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
%> @param   fieldType (string): select 'Total field', 'Scattered field'
%>          or 'Initial field'
%>
%> @param   scale (scalar): optional scale parameter for quiver plot. See
%           quiver documentation
%===============================================================================
function S = plot_poynting(ax,simulation,fieldType,scale)

if ~exist('scale','var')
    scale = 3;
end

hold(ax,'on')

pArr = simulation.input.particles.positionArray;
rArr = simulation.input.particles.radiusArray;
dims = simulation.output.fieldPointsArrayDims;

switch lower(fieldType)
    case 'total field'
        E = simulation.output.totalEField;
        H = simulation.output.totalHField;
    case 'scattered field'
        E = simulation.output.scatteredEField + simulation.output.internalEField;
        H = simulation.output.scatteredHField + simulation.output.internalHField;
    case 'initial field'
        E = simulation.output.initialEField;
        H = simulation.output.initialHField;
end

S = 0.5*real(cross(reshape(gather(E), [dims,3]), conj(reshape(gather(H), [dims,3]))));

fldPnts = reshape([simulation.output.fieldPoints(:,1), ...
                   simulation.output.fieldPoints(:,2), ...
                   simulation.output.fieldPoints(:,3)], [dims,3]);

if all(fldPnts(:,:,1) == fldPnts(1,1,1))        % fldPoints are on the yz plane
    perpdim = 1;                                % 1->x is the perp. direction
elseif all(fldPnts(:,:,2) == fldPnts(1,1,2))    % fldPoints are on the xz plane
    perpdim = 2;                                % 2->y is the perp. direction
elseif all(fldPnts(:,:,3) == fldPnts(1,1,3))    % fldPoints are on the xy plane
    perpdim = 3;                                % 3->z is the perp. direction
else
    error('fieldPoint must define an xy, xz or yz-like plane')
end

xy = setdiff([1,2,3], perpdim);                 % here xy are the in-plane dimensions
x = fldPnts(:,:,xy(1));
y = fldPnts(:,:,xy(2));

dist = abs(pArr(:,perpdim) - fldPnts(1,1,perpdim));% particle distances from xy plane
idx = find(dist<rArr);                          % find particles intersecting the plane
rArr(idx) = sqrt(rArr(idx).^2 - dist(idx).^2);  % overwrite radius of the intersection circle

quiver(x,y,S(:,:,xy(1)),S(:,:,xy(2)),scale), axis equal tight

for i=1:length(idx)
    rectangle(ax, ...
             'Position', [pArr(idx(i),xy)-rArr(idx(i)), [2,2]*rArr(idx(i))], ...
             'Curvature', [1 1], ...
             'FaceColor', 'none', ...
             'EdgeColor', [0,0,0], ...
             'LineWidth', 0.75)
end

labels = ['x'; 'y'; 'z'];
xlabel(labels(xy(1)))
ylabel(labels(xy(2)))

ax.DataAspectRatio = [1,1,1];
title([fieldType,', Poynting vector'])

hold(ax,'off')

end
