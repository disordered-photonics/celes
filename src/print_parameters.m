function print_parameters(simulation)

fprintf(1,'initial field: ');
if (simulation.input.initialField.beamWidth==0 || simulation.input.initialField.beamWidth==inf)
    fprintf(1,'plane wave.\n');
else
    fprintf(1,'Gaussian beam.\n');
end
fprintf(1,'wavelength: %.2f\n',simulation.input.wavelength);
fprintf(1,'number of particles: %i\n',simulation.input.particles.number);

