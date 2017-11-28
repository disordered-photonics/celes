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
%> @brief Generate a celes_simulation object from a model wizard app object
%>
%> @param       app (model_wizard application object)
%===============================================================================
function simulation = app2simulation(app)

simulation = celes_simulation;

% initialize fields of simulation object
simulation.input = celes_input;
simulation.numerics = celes_numerics;
simulation.tables = celes_tables;

% initialize fields of simulation.input object
simulation.input.wavelength = single(app.WavelengthEditField.Value);
simulation.input.mediumRefractiveIndex = single(app.BackgroundRefractiveIndexEditField.Value);
simulation.input.initialField = celes_initialField;
simulation.input.particles = celes_particles;

% initialize fields of simulation.input.particles object
simulation.input.particles.positionArray = single(app.particlePositions);
simulation.input.particles.refractiveIndexArray = app.particleRefractiveIndices;
simulation.input.particles.radiusArray = app.particleRadii;

% initialize fields of simulation.input.initialField object
simulation.input.initialField.amplitude = single(app.AmplitudeEditField.Value);
simulation.input.initialField.polarization = app.PolarizationDropDown.Value;
simulation.input.initialField.beamWidth = single(app.InitialBeamWaistEditField.Value);
simulation.input.initialField.focalPoint = single([app.InitialFocusXEditField.Value,app.InitialFocusYEditField.Value,app.InitialFocusZEditField.Value]);

% initialize fields of simulation.numerics object
simulation.numerics.lmax = app.lmaxEditField.Value;
simulation.numerics.particleDistanceResolution = single(app.LookupTableResolutionEditField.Value);
simulation.numerics.gpuFlag = app.ComputeOnGPUCheckBox.Value;
if app.CustomPolarGridCheckBox.Value
    simulation.numerics.polarAnglesArray = single(app.polarGrid);
else
    simulation.numerics.polarAnglesArray = single(0:app.PolarResolutionEditField.Value:180)*pi/180;
end
if app.CustomAzimuthalGridCheckBox.Value
    simulation.numerics.azimuthalAnglesArray = single(app.azimuthalGrid);
else
    simulation.numerics.azimuthalAnglesArray = single(0:app.AzimuthalResolutionEditField.Value:360)*pi/180;
end
simulation.numerics.solver = celes_solver;
simulation.numerics.solver.type = app.SolverTypeDropDown.Value;

% initialize fields of simulation.numerics.preconditioner object
simulation.numerics.solver.preconditioner = celes_preconditioner;
simulation.numerics.solver.preconditioner.type = app.PreconditionerDropDown.Value;
simulation.numerics.solver.preconditioner.partitionEdgeSizes = [app.PreconditionerPartitionXEditField.Value,app.PreconditionerPartitionYEditField.Value,app.PreconditionerPartitionZEditField.Value];

% initialize fields of simulation.numerics.solver object
simulation.numerics.solver.tolerance = app.SolverToleranceEditField.Value;
simulation.numerics.solver.maxIter = app.SolverMaxiterEditField.Value;

% initialize fields of simulation.output object
simulation.output = app2output(app);

end
