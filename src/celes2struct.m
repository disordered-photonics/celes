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

% ======================================================================
%> @brief Save relevant simulation fields to a struct
% ======================================================================

function s = celes2struct(obj)
validateattributes(obj, {'celes_simulation'}, {'nonempty'})

s.wavelength = obj.input.wavelength;
s.mediumRefractiveIndex = obj.input.mediumRefractiveIndex;
s.positionArray = obj.input.particles.positionArray;
s.refractiveIndexArray = obj.input.particles.refractiveIndexArray;
s.radiusArray = obj.input.particles.radiusArray;
s.polarAngle = obj.input.initialField.polarAngle;
s.azimuthalAngle = obj.input.initialField.azimuthalAngle;
s.polarization = obj.input.initialField.polarization;
s.beamWidth = obj.input.initialField.beamWidth;
s.focalPoint = obj.input.initialField.focalPoint;
s.lmax = obj.numerics.lmax;
s.polarAnglesArray = obj.numerics.polarAnglesArray;
s.azimuthalAnglesArray = obj.numerics.azimuthalAnglesArray;
s.gpuFlag = obj.numerics.gpuFlag;
s.particleDistanceResolution = obj.numerics.particleDistanceResolution;
s.solvertype = obj.numerics.solver.type;
s.tolerance = obj.numerics.solver.tolerance;
s.maxIter = obj.numerics.solver.maxIter;
s.restart = obj.numerics.solver.restart;
s.precondtype = obj.numerics.solver.preconditioner.type;
s.partitionEdgeSizes = obj.numerics.solver.preconditioner.partitionEdgeSizes;
end

