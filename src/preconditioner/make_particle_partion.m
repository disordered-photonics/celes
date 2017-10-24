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
%> @brief Generate a partition of the particles for the use by a
%>        blockdiagonal preconditioner
%>
%> @param   positionArray (Nx3 float array): particle positions [x,y,z]
%>
%> @param   edgeSizes (1x3 float array): edge sizes of the cuboids that
%>          define the partition
%>
%> @retval  partitioning (cell array): each cell of partitioning
%>          contains an array of inidices corresponding to spheres that
%>          fall into the same partition cuboid
%===============================================================================
function partitioning = make_particle_partion(positionArray,edgeSizes)

xarray = [(min(positionArray(:,1))-1):edgeSizes(1):(max(positionArray(:,1))+1),(max(positionArray(:,1))+1)];
yarray = [(min(positionArray(:,2))-1):edgeSizes(2):(max(positionArray(:,2))+1),(max(positionArray(:,2))+1)];
zarray = [(min(positionArray(:,3))-1):edgeSizes(3):(max(positionArray(:,3))+1),(max(positionArray(:,3))+1)];

[Nx,~,binx] = histcounts(positionArray(:,1),xarray);
[Ny,~,biny] = histcounts(positionArray(:,2),yarray);
[Nz,~,binz] = histcounts(positionArray(:,3),zarray);

partitioning = {};

for jx = 1:length(Nx)
    inx = find(binx == jx);
    for jy = 1:length(Ny)
        iny = find(biny == jy);
        inxy = intersect(inx, iny);
        for jz = 1:length(Nz)
            inz = find(binz == jz);
            inxyz = intersect(inxy, inz);
            if ~isempty(inxyz)
                partitioning{end+1} = inxyz;
            end
        end
    end 
end
end
