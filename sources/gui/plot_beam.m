%======================================================================
%> @brief Visualize beam width and focus
%>
%> @param       ax (axes object)
%> @param       app (model_wizard application object)
%======================================================================
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