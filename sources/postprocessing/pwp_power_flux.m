%======================================================================
%> @brief Compute the power flux for a field specified by a certain plane 
%> wave pattern
%>
%> The power flux is evaluated according to the formula
%> \f$ P = \frac{2\pi^2}{\omega k \mu_0}\sum_{j=1}^2
%> \int\dd{\alpha}\int\dd{\beta}\sin(\beta) \abs{g_j(\alpha,\beta)}^2 \f$,
%> see the \ref theory section.
%>
%> @param       pwp (celes_planeWavePattern): plane wave pattern of the
%>              field for which the power flux shall be evaluated
%>
%> @param       simulation (celes_simulation object): simulation object
%>              that stores all parameters of the simulation
%>
%> @param       Optional: direction (string): select 'forward' or
%>              'backward'. If specified, only partial waves propagating
%>              in the forward or backward z-direction are considered
%>
%> @retval      P (float): power flux
%======================================================================
function P=pwp_power_flux(pwp,simulation,varargin)


if isempty(varargin)
    beta = pwp.polarAngles;
    alpha = pwp.azimuthalAngles;
    g = pwp.expansionCoefficients;
else
    switch varargin{1}
        case 'forward'
            forward_idcs = (cos(pwp.polarAngles)>=0);
            beta = pwp.polarAngles(forward_idcs);
            alpha = pwp.azimuthalAngles;
            g = pwp.expansionCoefficients(:,forward_idcs);
        case 'backward'
            backward_idcs = (cos(pwp.polarAngles)<=0);
            beta = pwp.polarAngles(backward_idcs);
            alpha = pwp.azimuthalAngles;
            g = pwp.expansionCoefficients(:,backward_idcs);
    end
end

bintgrnd = trapz(alpha,abs(g).^2);
intgrl = trapz(beta, bintgrnd.*sin(beta));
P = 2*pi^2/(simulation.input.omega*simulation.input.k_medium)*real(intgrl);

