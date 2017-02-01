%======================================================================
%> @brief Show a 2D polar plot of the far field intensity
%>
%> @param       simulation (celes_simulation): simulation object
%>              containing the solution to plot
%>
%> @param       direction (string): select 'Forward intensity' or 
%>              'Backward intensity'
%>
%> @param       polarization (string): select 'TE', 'TM' or 'TE+TM'
%>
%> @param       fieldType (string): select 'Total field', 'Scattered field'
%>              or 'Initial field'
%======================================================================
function plot_intensity(simulation,direction,polarization,fieldType)
                    
switch fieldType
    case 'Total field'
        pwp = simulation.output.totalFieldPlaneWavePattern;
    case 'Scattered field'
        pwp = simulation.output.scatteredFieldPlaneWavePattern;
    case 'Initial field'
        pwp = initial_field_plane_wave_pattern(simulation);
end

switch polarization
    case 'TE'
        g2 = gather(abs(pwp{1}.expansionCoefficients).^2);
    case 'TM'
        g2 = gather(abs(pwp{2}.expansionCoefficients).^2);
    case 'TE+TM'
        g2 = gather(abs(pwp{1}.expansionCoefficients).^2+abs(pwp{2}.expansionCoefficients).^2);
end

switch direction
    case 'Forward intensity'
        forward_idcs = (cos(pwp{1}.polarAngles)>=0);
        g2 = g2(:,forward_idcs);
    case 'Backward intensity'
        backward_idcs = (cos(pwp{1}.polarAngles)<=0);
        g2 = g2(:,backward_idcs);
        g2 = g2(:,end:-1:1);
end

polarplot3d(double(g2).');
view([0,90])
set(gca,'DataAspectRatio',[1,1,1])