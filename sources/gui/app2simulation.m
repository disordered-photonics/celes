%======================================================================
%> @brief Generate a celes_simulation object from a model wizard app
%> object
%>
%> @param       app (model_wizard application object)
%======================================================================
function simulation = app2simulation(app)

simulation = celes_simulation;

% initialize fields of simulation object
simulation.input = celes_input;
simulation.numerics = celes_numerics;
simulation.tables = celes_tables;
simulation.output = celes_output;

% initialize fields of simulation.input object
simulation.input.wavelength = single(app.WavelengthEditField.Value);
simulation.input.mediumRefractiveIndex = single(app.BackgroundRefractiveIndexEditField.Value);
simulation.input.initialField = celes_initialField;
simulation.input.particles = celes_particles;

% initialize fields of simulation.input.particles object
simulation.input.particles.positionArray = single(app.particlePositions);
simulation.input.particles.refractiveIndex = single(app.ParticleRefractiveIndexEditField.Value)+1i*single(app.ParticleExtinctionCoefficientEditField.Value);
simulation.input.particles.radius = single(app.ParticleRadiusEditField.Value);

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
simulation.numerics.solver.monitor = app.SolverMonitorCheckBox.Value;

% initialize fields of simulation.output object
if app.FieldEvaluationCheckBox.Value
    dim1arr = app.FieldDim1MinEditField.Value:app.FieldDim1StepEditField.Value:app.FieldDim1MaxEditField.Value;
    dim2arr = app.FieldDim2MinEditField.Value:app.FieldDim2StepEditField.Value:app.FieldDim2MaxEditField.Value;
    [dim1,dim2]=meshgrid(dim1arr,dim2arr);
    dim3=dim1-dim1+app.FieldPlanePosEditField.Value;
    sz=size(dim1);
    switch app.FieldPlaneDropDown.Value
        case 'xz-Plane'
            simulation.output.fieldPoints=[dim1(:),dim3(:),dim2(:)];
        case 'yz-Plane'
            simulation.output.fieldPoints=[dim3(:),dim1(:),dim2(:)];
        case 'xy-Plane'
            simulation.output.fieldPoints=[dim1(:),dim2(:),dim3(:)];
    end
    simulation.output.fieldPointsArrayDims = sz;
end