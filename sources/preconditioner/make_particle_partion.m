%======================================================================
%> @brief Generate a partition of the particles for the use by a
%> blockdiagonal preconditioner
%>
%> @param       positionArray (Nx3 float array): particle positions [x,y,z]
%>
%> @param       edgeSizes (1x3 float array): edge sizes of the cuboids that
%>              define the partition
%>
%> @retval      partitioning (cell array): each cell of partitioning
%>              contains an array of inidices corresponding to spheres that
%>              fall into the same partition cuboid
%======================================================================
function partitioning = make_particle_partion(positionArray,edgeSizes)

xarray = [(min(positionArray(:,1))-1):edgeSizes(1):(max(positionArray(:,1))+1),(max(positionArray(:,1))+1)];
yarray = [(min(positionArray(:,2))-1):edgeSizes(2):(max(positionArray(:,2))+1),(max(positionArray(:,2))+1)];
zarray = [(min(positionArray(:,3))-1):edgeSizes(3):(max(positionArray(:,3))+1),(max(positionArray(:,3))+1)];
[Nx,~,binx] = histcounts(positionArray(:,1),xarray);
[Ny,~,biny] = histcounts(positionArray(:,2),yarray);
[Nz,~,binz] = histcounts(positionArray(:,3),zarray);
partitioning = {};
for jx=1:length(Nx)
    inx=find(binx==jx);
    for jy=1:length(Ny)
        iny=find(biny==jy);
        inxy=intersect(inx,iny);
        for jz=1:length(Nz)
            inz=find(binz==jz);
            inxyz=intersect(inxy,inz);
            if ~isempty(inxyz)
                partitioning{end+1}=inxyz;
            end
        end
    end
    
end