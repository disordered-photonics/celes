%======================================================================
%> @brief Read simulation parameters from celes_simulation object and
%> store them in model_wizard app object
%>
%> @param       simulation (celes_simulation)
%> @param       app (model_wizard application object)
%======================================================================
function simulation2app(simulation,app)

app.simulation = simulation;

% fields of simulation.input object
app.WavelengthEditField.Value = simulation.input.wavelength;
app.BackgroundRefractiveIndexEditField.Value = simulation.input.mediumRefractiveIndex;

% fields of simulation.input.particles object
app.particlePositions = simulation.input.particles.positionArray;
app.ParticleRefractiveIndexEditField.Value = real(simulation.input.particles.refractiveIndex);
app.ParticleExtinctionCoefficientEditField.Value = imag(simulation.input.particles.refractiveIndex);
app.ParticleRadiusEditField.Value = simulation.input.particles.radius;

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
app.SolverMonitorCheckBox.Value = simulation.numerics.solver.monitor;

app.FieldEvaluationCheckBox.Value = true;
app.FieldPlaneDropDown.Value = 'Custom';
app.fieldPoints = simulation.output.fieldPoints;