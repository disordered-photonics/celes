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
%> @brief Add circles representing the sphere positions to some axes object
%>
%> @param       ax (axes object): axes to plot in
%>
%> @param       positionArray (Nx3 float array): particle positions [x,y,z]
%>
%> @param       radiusArray (1xN float array): particle radius array
%>
%> @param       view (string): select 'view xy','view yz' or 'view xz'
%>
%> @param       fieldType (string): select 'Total field', 'Scattered field'
%>              or 'Initial field'
%>
%> @param       Optional: plotDepthInterval (2x1 float array): from where to where
%>              include spheres in the plot? Example: if the field points
%>              are all located in the xy-plane, potDepthInterval=[-200,200]
%>              will plot all spheres with a center z-coordinate in that
%>              interval
%======================================================================
function plot_spheres(ax,positionArray,radiusArray,refractiveIndexArray,view,varargin)

% normalize to 1
refractiveIndexArray = abs(refractiveIndexArray)/max(abs(refractiveIndexArray));

if isempty(varargin)
    plotDepthInterval = [-Inf,Inf];
else
    plotDepthInterval = varargin{1};
end

hold(ax,'on')

switch view
    case 'view xy'
        [~,idx]=sort(positionArray(:,3));
        for jS=1:length(positionArray(:,1))
            if positionArray(idx(jS),3)>plotDepthInterval(1) && positionArray(idx(jS),3)<plotDepthInterval(2)
                rectangle(ax,'Position',[positionArray(idx(jS),1:2)-[1,1]*radiusArray(jS),[2,2]*radiusArray(jS)],'Curvature',[1 1],'FaceColor',[0.5,0.5,refractiveIndexArray(jS)])
            end
        end
        ax.XLim = [min(positionArray(:,1))-3*max(radiusArray) , max(positionArray(:,1))+3*max(radiusArray)];
        ax.YLim = [min(positionArray(:,2))-3*max(radiusArray) , max(positionArray(:,2))+3*max(radiusArray)];
    case 'view xz'
        [~,idx]=sort(positionArray(:,2));
        for jS=1:length(positionArray(:,1))
            if positionArray(idx(jS),2)>plotDepthInterval(1) && positionArray(idx(jS),2)<plotDepthInterval(2)
                rectangle(ax,'Position',[positionArray(idx(jS),[1,3])-[1,1]*radiusArray(jS),[2,2]*radiusArray(jS)],'Curvature',[1 1],'FaceColor',[0.5,0.5,refractiveIndexArray(jS)])
            end
        end
        ax.XLim = [min(positionArray(:,1))-3*max(radiusArray) , max(positionArray(:,1))+3*max(radiusArray)];
        ax.YLim = [min(positionArray(:,3))-3*max(radiusArray) , max(positionArray(:,3))+3*max(radiusArray)];
    case 'view yz'
        [~,idx]=sort(positionArray(:,1));
        for jS=1:length(positionArray(:,1))
            if positionArray(idx(jS),1)>plotDepthInterval(1) && positionArray(idx(jS),1)<plotDepthInterval(2)
                rectangle(ax,'Position',[positionArray(idx(jS),2:3)-[1,1]*radiusArray(jS),[2,2]*radiusArray(jS)],'Curvature',[1 1],'FaceColor',[0.5,0.5,refractiveIndexArray(jS)])
            end
        end
        ax.XLim = [min(positionArray(:,2))-3*max(radiusArray) , max(positionArray(:,2))+3*max(radiusArray)];
        ax.YLim = [min(positionArray(:,3))-3*max(radiusArray) , max(positionArray(:,3))+3*max(radiusArray)];
end
ax.DataAspectRatio=[1,1,1];

hold(ax,'off')