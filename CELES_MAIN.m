% ======================================================================
%> @brief Run this script to start the simulation
% ======================================================================

% -------------------------------------------------------------------------
% do not edit
% -------------------------------------------------------------------------

% add all folders to the matlab search path
addpath(genpath('./src'))

% initialize the celes class instances
simulation = celes_simulation;
particles = celes_particles;
initialField = celes_initialField;
input = celes_input;
numerics = celes_numerics;
solver = celes_solver;
preconditioner = celes_preconditioner;
output = celes_output;


% -------------------------------------------------------------------------
% begin of user editable section - specify the simulation parameters here
% -------------------------------------------------------------------------

% import example sphere positions, radii and refractive indices
sphere_data = dlmread('examples/sphere_parameters.txt');

% radii of particles
% must be an array with one element per sphere
particles.radiusArray = sphere_data(:, 4);

% complex refractive index of particles, n+ik
% refractiveIndexArray must be same dimension as radiusArray
particles.refractiveIndexArray = sphere_data(:,5) + 1i*sphere_data(:, 6);

% positions of particles (in three-column format: x,y,z)
particles.positionArray = sphere_data(:,1:3);

% polar angle of incoming beam/wave, in radians (for Gaussian beams, 
% only 0 and pi are currently possible)
initialField.polarAngle = 0;

% azimuthal angle of incoming beam/wave, in radians
initialField.azimuthalAngle = 0;

% polarization of incoming beam/wave ('TE' or 'TM')
initialField.polarization = 'TE';

% width of beam waist (use 0 or inf for plane wave)
initialField.beamWidth = 2000;

% focal point 
initialField.focalPoint = [0,0,0];

% vacuum wavelength (same unit as particle positions and radius)
input.wavelength = 550;

% complex refractive index of surrounding medium
input.mediumRefractiveIndex = 1;

% maximal expansion order of scattered fields (around particle center)
numerics.lmax = 3;

% resolution of lookup table for spherical Hankel function (same unit as
% wavelength)
numerics.particleDistanceResolution = 1;

% use GPU for various calculations (deactivate if you experience GPU memory 
% problems - translation operator always runs on gpu, even if false)
numerics.gpuFlag = true;

% sampling of polar angles in the plane wave patterns (radians array)
numerics.polarAnglesArray = 0:pi/1e3:pi;

% sampling of azimuthal angles in the plane wave patterns (radians array)
numerics.azimuthalAnglesArray = 0:pi/1e3:2*pi;

% specify solver type (currently 'BiCGStab' or 'GMRES')
solver.type = 'GMRES';

% relative accuracy (solver stops when achieved)
solver.tolerance=5e-4;

% maximal number of iterations (solver stops if exceeded)
solver.maxIter=1000;

% restart parameter (only for GMRES)
solver.restart=1000;

% type of preconditioner (currently only 'blockdiagonal' and 'none'
% possible)
preconditioner.type = 'blockdiagonal';

% for blockdiagonal preconditioner: edge size of partitioning cuboids
preconditioner.partitionEdgeSizes = [1200,1200,1200];

% specify the points where to evaluate the electric near field (3-column
% array x,y,z)
[x,z] = meshgrid(-4000:50:4000,-3000:50:5000); y=x-x;
output.fieldPoints = [x(:),y(:),z(:)];

% dimensions of the array used to restore the original shape of x,y,z 
% for display of 2d-images
output.fieldPointsArrayDims = size(x);

% -------------------------------------------------------------------------
% end of user editable section
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% do not edit
% -------------------------------------------------------------------------

% combine all parameters to one simulation object
input.particles = particles;
input.initialField = initialField;
solver.preconditioner = preconditioner; 
numerics.solver = solver;
simulation.input = input;
simulation.numerics = numerics;
simulation.tables = celes_tables;
simulation.output = output;

% compute
simulation=simulation.run;
simulation=simulation.evaluatePower;
simulation=simulation.evaluateFields;

% display the particle positions
figure
plot_spheres(gca,simulation.input.particles.positionArray,simulation.input.particles.radiusArray,simulation.input.particles.refractiveIndexArray,'view xy')

% view near field
figure
plot_field(gca,simulation,'abs E','Total field',simulation.input.particles.radiusArray)
colormap(jet)
colorbar
caxis([0,1.2])

fprintf('transmitted power: %f %%\n',simulation.output.totalFieldForwardPower/simulation.output.initialFieldPower*100)
fprintf('reflected power: %f %%\n',simulation.output.totalFieldBackwardPower/simulation.output.initialFieldPower*100)
