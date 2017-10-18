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
%===============================================================================
function plot_intensity(simulation,direction,polarization,fieldType)

switch lower(fieldType)
    case 'total field'
        pwp = simulation.output.totalFieldPlaneWavePattern;
    case 'scattered field'
        pwp = simulation.output.scatteredFieldPlaneWavePattern;
    case 'initial field'
        pwp = initial_field_plane_wave_pattern(simulation);
end

switch lower(polarization)
    case 'te'
        g2 = gather(abs(pwp{1}.expansionCoefficients).^2);
    case 'tm'
        g2 = gather(abs(pwp{2}.expansionCoefficients).^2);
    case 'te+tm'
        g2 = gather(abs(pwp{1}.expansionCoefficients).^2+ ...
                    abs(pwp{2}.expansionCoefficients).^2);
end

switch lower(direction)
    case 'forward intensity'
        forward_idcs = (cos(pwp{1}.polarAngles)>=0);
        g2 = g2(:,forward_idcs);
    case 'backward intensity'
        backward_idcs = (cos(pwp{1}.polarAngles)<=0);
        g2 = g2(:,backward_idcs);
        g2 = g2(:,end:-1:1);
end

polarplot3d(double(g2).');
view([0,90])
set(gca,'DataAspectRatio',[1,1,1])
end
