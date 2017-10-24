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
%> @brief Read simulation parameters from celes_simulation object and
%>        store them in model_wizard app object
%>
%> @param       simulation (celes_simulation)
%> @param       app (model_wizard application object)
%===============================================================================
function simulation2app(simulation,app)

app.simulation = simulation;

% fields of simulation.input object
app.WavelengthEditField.Value = simulation.input.wavelength;
app.BackgroundRefractiveIndexEditField.Value = simulation.input.mediumRefractiveIndex;

% fields of simulation.input.particles object
app.particlePositions = simulation.input.particles.positionArray;
app.particleRefractiveIndices = simulation.input.particles.refractiveIndexArray;
app.particleRadii = simulation.input.particles.radiusArray;

% fields of simulation.input.initialField object
app.AmplitudeEditField.Value = simulation.input.initialField.amplitude;
app.PolarizationDropDown.Value = simulation.input.initialField.polarization;
app.InitialBeamWaistEditField.Value = simulation.input.initialField.beamWidth;
foc = simulation.input.initialField.focalPoint;
app.InitialFocusXEditField.Value = foc(1);
app.InitialFocusYEditField.Value = foc(2);
app.InitialFocusZEditField.Value = foc(3);

% fields of simulation.numerics object
app.lmaxEditField.Value = simulation.numerics.lmax;
app.LookupTableResolutionEditField.Value = simulation.numerics.particleDistanceResolution;
app.ComputeOnGPUCheckBox.Value = simulation.numerics.gpuFlag;
app.CustomPolarGridCheckBox.Value = true;
app.PolarResolutionEditField.Enable = 'off';
app.polarGrid = simulation.numerics.polarAnglesArray;
app.CustomAzimuthalGridCheckBox.Value = true;
app.AzimuthalResolutionEditField.Enable = 'off';
app.azimuthalGrid = simulation.numerics.azimuthalAnglesArray;

app.SolverTypeDropDown.Value = simulation.numerics.solver.type;
app.PreconditionerDropDown.Value = simulation.numerics.solver.preconditioner.type;
edgsiz = simulation.numerics.solver.preconditioner.partitionEdgeSizes;
app.PreconditionerPartitionXEditField.Value = edgsiz(1);
app.PreconditionerPartitionYEditField.Value = edgsiz(2);
app.PreconditionerPartitionZEditField.Value = edgsiz(3);

app.SolverToleranceEditField.Value = simulation.numerics.solver.tolerance;
app.SolverMaxiterEditField.Value = simulation.numerics.solver.maxIter;

app.FieldEvaluationCheckBox.Value = false;
end
