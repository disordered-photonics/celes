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

%======================================================================
%> @brief Show an image-plot of the near field
%>
%> The sphere positions are displayed as circles.
%>
%> @param       ax (axes object): axes to plot in
%>
%> @param       simulation (celes_simulation): simulation object
%>              containing the solution to plot
%>
%> @param       component (string): select 'real Ex', 'real Ey', 'real Ez'
%>              or 'abs E'
%>
%> @param       fieldType (string): select 'Total field', 'Scattered field'
%>              or 'Initial field'
%>
%> @param       radiusArray (1xN float array): radius array of the spheres
%>
%> @param       Optional: plotDepthInterval (2x1 float array): from where to where
%>              include spheres in the plot? Example: if the field points
%>              are all located in the xy-plane, potDepthInterval=[-200,200]
%>              will plot all spheres with a center z-coordinate in that
%>              interval
%======================================================================
function plot_field(ax,simulation,component,fieldType,radiusArray,varargin)

hold(ax,'on')

switch fieldType
    case 'Total field'
        E = simulation.output.totalField;
    case 'Scattered field'
        E = simulation.output.scatteredField + simulation.output.internalField;
    case 'Initial field'
        E = simulation.output.initialField;
end

x = simulation.output.fieldPoints(:,1);
y = simulation.output.fieldPoints(:,2);
z = simulation.output.fieldPoints(:,3);

if all(x==x(1))
    view = 'yz';
    y = reshape(y,simulation.output.fieldPointsArrayDims);
    z = reshape(z,simulation.output.fieldPointsArrayDims);
elseif all(y==y(1))
    view = 'xz';
    x = reshape(x,simulation.output.fieldPointsArrayDims);
    z = reshape(z,simulation.output.fieldPointsArrayDims);
elseif all(z==z(1))
    view = 'xy';
    x = reshape(x,simulation.output.fieldPointsArrayDims);
    y = reshape(y,simulation.output.fieldPointsArrayDims);
end   

switch component
    case 'real Ex'
        fld = reshape(gather(real(E(:,1))),simulation.output.fieldPointsArrayDims);
    case 'real Ey'
        fld = reshape(gather(real(E(:,2))),simulation.output.fieldPointsArrayDims);
    case 'real Ez'
        fld = reshape(gather(real(E(:,3))),simulation.output.fieldPointsArrayDims);
    case 'abs E'
        fld = reshape(gather(sqrt(abs(E(:,1)).^2 + abs(E(:,2)).^2 + abs(E(:,3)).^2)),simulation.output.fieldPointsArrayDims);
end

if isempty(varargin)
    plotDepthInterval = [-Inf,Inf];
else
    plotDepthInterval = varargin{1};
end

hold(ax,'on')

positionArray = simulation.input.particles.positionArray;

switch view
    case 'xy'
        imagesc(x(1,:),y(:,1),fld)
        xlabel('x')
        ylabel('y')
        [~,idx]=sort(positionArray(:,3));
        for jS=1:length(positionArray(:,1))
            if positionArray(idx(jS),3)-z(1)>plotDepthInterval(1) && positionArray(idx(jS),3)-z(1)<plotDepthInterval(2)
                rectangle(ax,'Position',[positionArray(idx(jS),1:2)-[1,1]*radiusArray(jS),[2,2]*radiusArray(jS)],'Curvature',[1 1],'FaceColor','none','EdgeColor',[1,1,1])
            end
        end
    case 'xz'
        imagesc(x(1,:),z(:,1),fld)
        xlabel('x')
        ylabel('z')
        [~,idx]=sort(positionArray(:,2));
        for jS=1:length(positionArray(:,1))
            if positionArray(idx(jS),2)>plotDepthInterval(1) && positionArray(idx(jS),2)<plotDepthInterval(2)
                rectangle(ax,'Position',[positionArray(idx(jS),[1,3])-[1,1]*radiusArray(jS),[2,2]*radiusArray(jS)],'Curvature',[1 1],'FaceColor','none','EdgeColor',[1,1,1])
            end
        end
    case 'yz'
        imagesc(y(1,:),z(:,1),fld)
        xlabel('y')
        ylabel('z')
        [~,idx]=sort(positionArray(:,1));
        for jS=1:length(positionArray(:,1))
            if positionArray(idx(jS),1)>plotDepthInterval(1) && positionArray(idx(jS),1)<plotDepthInterval(2)
                rectangle(ax,'Position',[positionArray(idx(jS),2:3)-[1,1]*radiusArray(jS),[2,2]*radiusArray(jS)],'Curvature',[1 1],'FaceColor','none','EdgeColor',[1,1,1])
            end
        end
end
ax.DataAspectRatio=[1,1,1];
axis tight
title([fieldType,', ',component])
hold(ax,'off')