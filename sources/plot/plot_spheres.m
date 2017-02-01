%======================================================================
%> @brief Add circles representing the sphere positions to some axes object
%>
%> @param       ax (axes object): axes to plot in
%>
%> @param       positionArray (Nx3 float array): particle positions [x,y,z]
%>
%> @param       radius (float): particle radius
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
function plot_spheres(ax,positionArray,radius,view,varargin)

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
                rectangle(ax,'Position',[positionArray(idx(jS),1:2)-[1,1]*radius,[2,2]*radius],'Curvature',[1 1],'FaceColor',[0.5,0.5,0.5])
            end
        end
        ax.XLim = [min(positionArray(:,1))-3*radius , max(positionArray(:,1))+3*radius];
        ax.YLim = [min(positionArray(:,2))-3*radius , max(positionArray(:,2))+3*radius];
    case 'view xz'
        [~,idx]=sort(positionArray(:,2));
        for jS=1:length(positionArray(:,1))
            if positionArray(idx(jS),2)>plotDepthInterval(1) && positionArray(idx(jS),2)<plotDepthInterval(2)
                rectangle(ax,'Position',[positionArray(idx(jS),[1,3])-[1,1]*radius,[2,2]*radius],'Curvature',[1 1],'FaceColor',[0.5,0.5,0.5])
            end
        end
        ax.XLim = [min(positionArray(:,1))-3*radius , max(positionArray(:,1))+3*radius];
        ax.YLim = [min(positionArray(:,3))-3*radius , max(positionArray(:,3))+3*radius];
    case 'view yz'
        [~,idx]=sort(positionArray(:,1));
        for jS=1:length(positionArray(:,1))
            if positionArray(idx(jS),1)>plotDepthInterval(1) && positionArray(idx(jS),1)<plotDepthInterval(2)
                rectangle(ax,'Position',[positionArray(idx(jS),2:3)-[1,1]*radius,[2,2]*radius],'Curvature',[1 1],'FaceColor',[0.5,0.5,0.5])
            end
        end
        ax.XLim = [min(positionArray(:,2))-3*radius , max(positionArray(:,2))+3*radius];
        ax.YLim = [min(positionArray(:,3))-3*radius , max(positionArray(:,3))+3*radius];
end
ax.DataAspectRatio=[1,1,1];

hold(ax,'off')