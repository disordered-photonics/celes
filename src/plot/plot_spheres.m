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
%> @brief Plot spheres in 3D, spheres' colors depend on their refractive index
%>
%> @param   ax (axes object): axes to plot in
%>
%> @param   positionArray (Nx3 float array): particle positions [x,y,z]
%>
%> @param   radiusArray (1xN float array): particle radius array
%>
%> @param   refractiveIndexArray (1xN float array): particle refractive index array
%>
%> @param   sphereResolution (scalar int): optional resolution for sphere()
%===============================================================================
function plot_spheres(ax, positionArray, radiusArray, refractiveIndexArray, sphereResolution)

if ~exist('sphereResolution','var')
    sphres = 16; % looks sufficiently like a sphere
else
    sphres = sphereResolution;
end
[vx, vy, vz] = sphere(sphres); % sphere base model

% normalize to 1
realRInorm = real(refractiveIndexArray)/max(abs((refractiveIndexArray)));
imagRInorm = imag(refractiveIndexArray)/max(abs(imag(refractiveIndexArray)));
if isnan(imagRInorm)
    imagRInorm = imag(refractiveIndexArray);
end

realRInorm = repmat(repelem(realRInorm, sphres+2),[1,sphres+1]);
imagRInorm = repmat(repelem(imagRInorm, sphres+2),[1,sphres+1]);

hold(ax,'on')

% plotting just one merged surface is incredibly faster
xmerged = zeros(size(positionArray,1)*(sphres+2), sphres+1);
ymerged = zeros(size(positionArray,1)*(sphres+2), sphres+1);
zmerged = zeros(size(positionArray,1)*(sphres+2), sphres+1);

for jS = 1:size(positionArray,1)
    idx = (1:sphres+2)+(jS-1)*(sphres+2);
    xmerged(idx,:) = [vx*radiusArray(jS) + positionArray(jS,1); NaN*ones(1,sphres+1)];
    ymerged(idx,:) = [vy*radiusArray(jS) + positionArray(jS,2); NaN*ones(1,sphres+1)];
    zmerged(idx,:) = [vz*radiusArray(jS) + positionArray(jS,3); NaN*ones(1,sphres+1)];
end

% color each sphere with a unique shade
CO(:,:,1) = imagRInorm; % red
CO(:,:,2) = zeros(size(zmerged)); % green
CO(:,:,3) = realRInorm; % blue

surf(xmerged, ymerged, zmerged, CO, 'LineStyle', 'none');

light; lighting gouraud
ax.DataAspectRatio = [1,1,1];
xlabel('x'); ylabel('y'); zlabel('z');

hold(ax,'off')
end
