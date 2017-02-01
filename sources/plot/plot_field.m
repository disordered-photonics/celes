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
%> @param       radius (float): radius of the spheres
%>
%> @param       Optional: plotDepthInterval (2x1 float array): from where to where
%>              include spheres in the plot? Example: if the field points
%>              are all located in the xy-plane, potDepthInterval=[-200,200]
%>              will plot all spheres with a center z-coordinate in that
%>              interval
%======================================================================
function plot_field(ax,simulation,component,fieldType,radius,varargin)

hold(ax,'on')

switch fieldType
    case 'Total field'
        E = simulation.output.totalField;
    case 'Scattered field'
        E = simulation.output.scatteredField;
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
                rectangle(ax,'Position',[positionArray(idx(jS),1:2)-[1,1]*radius,[2,2]*radius],'Curvature',[1 1],'FaceColor',[0.5,0.5,0.5])
            end
        end
    case 'xz'
        imagesc(x(1,:),z(:,1),fld)
        xlabel('x')
        ylabel('z')
        [~,idx]=sort(positionArray(:,2));
        for jS=1:length(positionArray(:,1))
            if positionArray(idx(jS),2)>plotDepthInterval(1) && positionArray(idx(jS),2)<plotDepthInterval(2)
                rectangle(ax,'Position',[positionArray(idx(jS),[1,3])-[1,1]*radius,[2,2]*radius],'Curvature',[1 1],'FaceColor',[0.5,0.5,0.5])
            end
        end
    case 'yz'
        imagesc(y(1,:),z(:,1),fld)
        xlabel('y')
        ylabel('z')
        [~,idx]=sort(positionArray(:,1));
        for jS=1:length(positionArray(:,1))
            if positionArray(idx(jS),1)>plotDepthInterval(1) && positionArray(idx(jS),1)<plotDepthInterval(2)
                rectangle(ax,'Position',[positionArray(idx(jS),2:3)-[1,1]*radius,[2,2]*radius],'Curvature',[1 1],'FaceColor',[0.5,0.5,0.5])
            end
        end
end
ax.DataAspectRatio=[1,1,1];
axis tight
title([fieldType,', ',component])
hold(ax,'off')