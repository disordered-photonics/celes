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
%> @brief Visualize beam width and focus
%>
%> @param       ax (axes object)
%> @param       app (model_wizard application object)
%===============================================================================
function plot_beam(ax,app)

beam_color = [1,0.8,0];
hold(ax,'on');

w0 = app.InitialBeamWaistEditField.Value;
wl = app.WavelengthEditField.Value;
xG = app.InitialFocusXEditField.Value;
yG = app.InitialFocusYEditField.Value;
zG = app.InitialFocusZEditField.Value;
switch app.ViewDropDown.Value
    case 'view xy'
        phi=0:0.01:2*pi;
        rho=w0;
        plot(ax,rho*cos(phi)+xG,rho*sin(phi)+yG,'LineWidth',2,'Color',beam_color)
    case 'view yz'
        z = ax.YLim(1):ax.YLim(2);
        y1 = w0*sqrt(1+((z-zG)*wl/pi/w0^2).^2)+yG;
        y2 = -w0*sqrt(1+((z-zG)*wl/pi/w0^2).^2)+yG;
        plot(ax,y1,z,'LineWidth',2,'Color',beam_color)
        plot(ax,y2,z,'LineWidth',2,'Color',beam_color)
    case 'view xz'
        z = ax.YLim(1):ax.YLim(2);
        x1 = w0*sqrt(1+((z-zG)*wl/pi/w0^2).^2)+xG;
        x2 = -w0*sqrt(1+((z-zG)*wl/pi/w0^2).^2)+xG;
        plot(ax,x1,z,'LineWidth',2,'Color',beam_color)
        plot(ax,x2,z,'LineWidth',2,'Color',beam_color)
end
hold(ax,'off')
end
