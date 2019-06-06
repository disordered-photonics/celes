% =========================================================================
%> @brief Run this script to start the simulation
% =========================================================================

%% add all folders to the MATLAB search path
addpath(genpath('./src'))

%% import example file with sphere positions, radii and refractive indices
data = dlmread('examples/sphere_parameters.txt');

%% initialize the CELES class instances

% initialize particle class instance
%   - positionArray:        Nx3 array (float) in [x,y,z] format
%   - refractiveIndexArray: Nx1 array (complex) of complex refractive indices
%   - radiusArray:          Nx1 array (float) of sphere radii
particles = celes_particles('positionArray',        data(:,1:3), ...
                            'refractiveIndexArray', data(:,5)+1i*data(:,6), ...
                            'radiusArray',          data(:,4) ...
                            );

% initialize initial field class instance
%   - polarAngle:           scalar (float) polar angle of incoming beam/wave,
%                           in radians. for Gaussian beams, only 0 or pi are
%                           currently possible
%   - azimuthalAngle:       scalar (float) azimuthal angle of incoming
%                           beam/wave, in radians
%   - polarization:         string (char) polarization of incoming beam/wave
%                           ('TE' or 'TM')
%   - beamWidth:            scalar (float) width of beam waist. use 0 or inf
%                           for plane wave
%   - focalPoint:           1x3 array (float) focal point
initialField = celes_initialField('polarAngle',     0, ...
                                  'azimuthalAngle', 0, ...
                                  'polarization',   'TE', ...
                                  'beamWidth',      2000, ...
                                  'focalPoint',     [0,0,0] ...
                                  );

% initialize input class instance
%   - wavelength:           scalar (float) vacuum wavelength, same unit as
%                           particle positions and radii
%   - mediumRefractiveIndex: scalar (complex) refractive index of environment
%   - particles:            valid instance of celes_particles class
%   - initialField:         valid instance of celes_initialField class
input = celes_input('wavelength',                   550, ...
                    'mediumRefractiveIndex',        1, ...
                    'particles',                    particles, ...
                    'initialField',                 initialField ...
                    );

% initialize preconditioner class instance
%   - type:                 string (char) type of preconditioner (currently
%                           only 'blockdiagonal' and 'none' available)
%   - partitionEdgeSizes    1x3 array (float) edge size of partitioning cuboids
%                           (applies to 'blockdiagonal' type only)
precnd = celes_preconditioner('type',               'blockdiagonal', ...
                              'partitionEdgeSizes', [1200,1200,1200] ...
                              );

% initialize solver class instance
%   - type:                 string (char) solver type (currently 'BiCGStab' or
%                           'GMRES' are implemented)
%   - tolerance:            scalar (float) target relative accuracy of solution
%   - maxIter:              scalar (int) maximum number of iterations allowed
%   - restart:              scalar (int) restart parameter (applies only to
%                           GMRES solver)
%   - preconditioner:       valid instance of celes_preconditioner class
solver = celes_solver('type',                       'GMRES', ...
                      'tolerance',                  5e-4, ...
                      'maxIter',                    1000, ...
                      'restart',                    200, ...
                      'preconditioner',             precnd);

% initialize numerics class instance
%   - lmax:                 scalar (int) maximal expansion order of scattered
%                           fields around particle center
%   - polarAnglesArray:     1xN array (float) sampling of polar angles in the
%                           plane wave patterns, in radians
%   - azimuthalAnglesArray: sampling of azimuthal angles in the plane wave
%                           patterns, in radians
%   - gpuFlag:              scalar (bool) set to false if you experience GPU
%                           memory problems at evaluation time (translation
%                           operator always runs on GPU, though)
%   - particleDistanceResolution: scalar (float) resolution of lookup table for
%                           spherical Hankel function (same unit as wavelength)
%   - solver:               valid instance of celes_solver class
numerics = celes_numerics('lmax',                   3, ...
                          'polarAnglesArray',       0:pi/5e3:pi, ...
                          'azimuthalAnglesArray',   0:pi/1e2:2*pi, ...
                          'gpuFlag',                true, ...
                          'particleDistanceResolution', 1, ...
                          'solver',                 solver);

% define a grid of points where the field will be evaulated
[x,z] = meshgrid(-4000:50:4000, -3000:50:5000); y = zeros(size(x));
% initialize output class instance
%   - fieldPoints:          Nx3 array (float) points where to evaluate the
%                           electric near field
%   - fieldPointsArrayDims: 1x2 array (int) dimensions of the array, needed to
%                           recompose the computed field as a n-by-m image
output = celes_output('fieldPoints',                [x(:),y(:),z(:)], ...
                      'fieldPointsArrayDims',       size(x));

% initialize simulation class instance
%   - input:                valid instance of celes_input class
%   - numerics:             valid instance of celes_input class
%   - output:               valid instance of celes_output class
simul = celes_simulation('input',                   input, ...
                         'numerics',                numerics, ...
                         'output',                  output);

%% run simulation
simul.run;

% evaluate transmitted and reflected power
simul.evaluatePower;
fprintf('transmitted power: \t%.4f %%\n', ...
    simul.output.totalFieldForwardPower/simul.output.initialFieldPower*100)
fprintf('reflected power: \t%.4f %%\n', ...
    simul.output.totalFieldBackwardPower/simul.output.initialFieldPower*100)

% evaluate field at output.fieldPoints
simul.evaluateFields;

%% plot results
% display particles
figure('Name','Particle positions','NumberTitle','off');
plot_spheres(gca,simul.input.particles.positionArray, ...
                 simul.input.particles.radiusArray, ...
                 simul.input.particles.refractiveIndexArray)

% plot near field
figure('Name','Near-field cross-cut','NumberTitle','off');
plot_field(gca,simul,'abs E','Total field')
caxis([0,2])

% % export animated gif
% figure('Name','Animated near-field cross-cut','NumberTitle','off');
% plot_field(gca,simul,'real Ey','Total field','Ey_total.gif')
