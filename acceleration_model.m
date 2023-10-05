clc
close all
clear all
format long

%% Global Constants

muo = 4*pi*10^(-7);
au = 1.66053892*10^(-27);
e = 1.60217657*10^(-19);
me = 9.10938291*10^(-31); 
g = 9.80665; % gravity

%% Engine Parameters (Per Engine)

% These are the parameters required by the user of the model
% Two main parameters we will control.

Id = 23000; % Current intensity [A]
mdot = 0.006; % Mass flow rate [kg/s]

% Geometry of the engine
Rc = 0.01; % Cathode radius [m]
Ra = 0.05; % Anode radius [m]
Lz = 0.1; % Channel length [m]

% Propellant
m_propellant = 39.948; % Propellant mass un au units
mi = m_propellant*au; % Actual mass of the propellant

% Approximations for the collisional frequency
Te = 5; % Temperature of the electrons [eV]

%% Geometry 

% This part of the code calculates the remaining geometrical parameters to
% know.

Ly = pi*(Ra+Rc); % 1D plan length [m]
Lr = Ra-Rc; % Radial plane length[m]
A = Ly*Lr; % Area [m^2]

%% PARAMETERS

% Calculates the relevant parameters of the engine using direct formulas of
% the axial model.

Bo = Id*muo/Ly; % Azimutal (y) magnetic field [T]
Go = mdot/A; % Mass flux [kg/s·m^2]
uE = Bo^2/(2*muo*Go); % Exit velocity [m/s]
Isp = uE/g; % Specific Impulse [s]
F = Id^2*muo*Lr/(2*Ly); % Thrust [N]
F_ln = muo*Id^2*log(Ra/Rc)/(4*pi); % Thrust calculated by Maecker's law [N]
nE = Go/(mi*uE); % Number of electrons
lnLAM_E = 9+0.5*(log(((10^18)/nE)*(Te)^3)); % Actual lnLAMBDA at the exit.
nue = (nE/(10^18))*(1/Te)^(3/2)*lnLAM_E*2.9*10^(6); % Collisional frequency
Q_ei = nue*2/nE; % Collisional rate
SmII = e^2*nE/(me*nue); % Parallel conductivity
wco = e*Bo/me; % Cyclotron frequency
Smpp = SmII*nue^2/(wco^2+nue^2); % Perpendicular conductivity
Rmo = SmII*muo*Lz*uE; % Characteristic magnetic reynolds number
uEBo = uE*Bo; % Induced E field
Puse = F*uE*0.5; % Power of the jet [W]

% Function that evaluates the engine performance at that reynolds
[Eopt,b,zn] = Bisection_method_E_nondim_b_u(Rmo);

%% Dimensionalisation of Model

B = b*Bo;
un = 1-b.^2;
u = un*uE;
uB = u.*B;
ub = un.*b;
z = zn*Lz;
s = size(z,1);
Eplot = linspace(Eopt,Eopt,s);
E = Eopt*uE*Bo;
Vd = E*Lr;
Ez = wco.*b.*Eplot'./nue;
j = SmII*(E-uB);
Eta_p = uE*Bo/(4*E);
PQN = Lr*E*Id;

%% Plotting

set(0,'DefaultAxesFontSize',22)
set(0,'DefaultAxesFontName','Vijaya')
set(0, 'DefaultAxesLineStyleOrder', '-');
set(0, 'DefaultAxesColorOrder', [0.0 0.0 0.0; 0.4 0.4 0.4; 0.6 0.6 0.6]);
plot(zn,Eplot','k',zn,b,'r',zn,un,'b',zn,ub,'m','LineWidth',2)
axis tight
set(gcf, 'Units', 'centimeters');
afFigurePosition = [10 6 18 12]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'Units','normalized','Position',[0.15 0.2 0.75 0.7]);
iFontSize = 20;
strFontUnit = 'points'; % [{points} | normalized | inches | centimeters | pixels]
strFontName = 'Times'; % [Times | Courier | ] TODO complete the list
strFontWeight = 'normal'; % [light | {normal} | demi | bold]
strFontAngle = 'normal'; % [{normal} | italic | oblique] ps: only for axes
strInterpreter = 'latex'; % [{tex} | latex]
fLineWidth = 1; % width of the line of the axes
set(gca, ...
... 'Position', [1 1 20 10], ... TODO
... 'OuterPosition', [1 1 20 10], ... TODO
...
'XGrid', 'off', ... [on | {off}]
'YGrid', 'off', ... [on | {off}]
'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
'XMinorGrid', 'off' , ... [on | {off}]
'YMinorGrid', 'off', ... [on | {off}]
'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
...
...'XTick', 0:0.1:100, ... ticks of x axis
...'YTick', 0:1:10, ... ticks of y axis
...'XTickLabel', {'-1','0','1'}, ...
...'YTickLabel', {'-1','0','1'}, ...
'XMinorTick', 'off', ... [on | {off}]
'YMinorTick', 'off', ... [on | {off}]
'TickDir', 'out', ... [{in} | out] inside or outside (for 2D)
'TickLength', [.01 .01], ... length of the ticks
...
'XColor', [.1 .1 .1], ... color of x axis
'YColor', [.1 .1 .1], ... color of y axis
'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
'XLim', [0 1], ... limits for the x-axis
'YLim', [0 1], ... limits for the y-axis
...
'FontName', strFontName, ... kind of fonts of labels
'FontSize', iFontSize, ... size of fonts of labels
'FontUnits', strFontUnit, ... units of the size of fonts
'FontWeight', strFontWeight, ... weight of fonts of labels
'FontAngle', strFontAngle, ... inclination of fonts of labels
...
'LineWidth', fLineWidth); % width of the line of the axes
% fonts properties
xlab = xlabel('$\bar{ z } [ - ]$');
set(xlab,'interpreter','Latex','FontSize',20)
leyenda= legend('$\bar{Ex}$','b','u','$u \times b$');
set(leyenda,'interpreter','Latex','FontSize',14,'LineWidth',2,'Location','North')
ylab = ylabel('$\bar{E}$ [ - ]');
set(ylab,'interpreter','Latex','FontSize',20)
figure
set(0,'DefaultAxesFontSize',12)
set(gcf, 'Units', 'centimeters');
afFigurePosition = [10 6 18 12]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'Units','normalized','Position',[0.15 0.2 0.75 0.7]);
iFontSize = 20;
strFontUnit = 'points'; % [{points} | normalized | inches | centimeters | pixels]
strFontName = 'Times'; % [Times | Courier | ] TODO complete the list

strFontWeight = 'normal'; % [light | {normal} | demi | bold]
strFontAngle = 'normal'; % [{normal} | italic | oblique] ps: only for axes
strInterpreter = 'latex'; % [{tex} | latex]
fLineWidth = 1; % width of the line of the axes
subplot(2,2,1)
plot(z,B,'LineWidth',2)
set(gca, ...
... 'Position', [1 1 20 10], ... TODO
... 'OuterPosition', [1 1 20 10], ... TODO
...
'XGrid', 'off', ... [on | {off}]
'YGrid', 'off', ... [on | {off}]
'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
'XMinorGrid', 'off' , ... [on | {off}]
'YMinorGrid', 'off', ... [on | {off}]
'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
...
...'XTick', 0:0.1:100, ... ticks of x axis
...'YTick', 0:1:10, ... ticks of y axis
...'XTickLabel', {'-1','0','1'}, ...
...'YTickLabel', {'-1','0','1'}, ...
'XMinorTick', 'off', ... [on | {off}]
'YMinorTick', 'off', ... [on | {off}]
'TickDir', 'out', ... [{in} | out] inside or outside (for 2D)
'TickLength', [.01 .01], ... length of the ticks
...
'XColor', [.1 .1 .1], ... color of x axis
'YColor', [.1 .1 .1], ... color of y axis
'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
...'XLim', [0 1], ... limits for the x-axis
...'YLim', [0 1], ... limits for the y-axis
...
'FontName', strFontName, ... kind of fonts of labels
'FontSize', iFontSize, ... size of fonts of labels
'FontUnits', strFontUnit, ... units of the size of fonts
'FontWeight', strFontWeight, ... weight of fonts of labels
'FontAngle', strFontAngle, ... inclination of fonts of labels
...
'LineWidth', fLineWidth); % width of the line of the axes
%fonts properties
axis tight
xlab = xlabel('z [ m ]');
set(xlab,'interpreter','Latex','FontSize',20)
ylab = ylabel('$B_{y}$ [ T ]');
set(ylab,'interpreter','Latex','FontSize',14)
axis tight
subplot(2,2,2)
plot(z,u,'LineWidth',2)
axis tight
set(gca, ...
... 'Position', [1 1 20 10], ... TODO
... 'OuterPosition', [1 1 20 10], ... TODO
...
'XGrid', 'off', ... [on | {off}]
'YGrid', 'off', ... [on | {off}]
'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
'XMinorGrid', 'off' , ... [on | {off}]
'YMinorGrid', 'off', ... [on | {off}]
'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
...
...'XTick', 0:0.1:100, ... ticks of x axis
...'YTick', 0:1:10, ... ticks of y axis
...'XTickLabel', {'-1','0','1'}, ...
...'YTickLabel', {'-1','0','1'}, ...
'XMinorTick', 'off', ... [on | {off}]
'YMinorTick', 'off', ... [on | {off}]
'TickDir', 'out', ... [{in} | out] inside or outside (for 2D)
'TickLength', [.01 .01], ... length of the ticks
...
'XColor', [.1 .1 .1], ... color of x axis
'YColor', [.1 .1 .1], ... color of y axis
'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
...'XLim', [0 1], ... limits for the x-axis
...'YLim', [0 1], ... limits for the y-axis
...
'FontName', strFontName, ... kind of fonts of labels
'FontSize', iFontSize, ... size of fonts of labels
'FontUnits', strFontUnit, ... units of the size of fonts
'FontWeight', strFontWeight, ... weight of fonts of labels
'FontAngle', strFontAngle, ... inclination of fonts of labels
...
'LineWidth', fLineWidth); % width of the line of the axes
% fonts properties
xlab = xlabel('z [ m ]');
set(xlab,'interpreter','Latex','FontSize',20)
ylab = ylabel('$u_{z} [ m s^{-1}]$');
set(ylab,'interpreter','Latex','FontSize',14)
subplot(2,2,3)
plot(z,uB,'LineWidth',2)
axis tight
set(gca, ...
... 'Position', [1 1 20 10], ... TODO
... 'OuterPosition', [1 1 20 10], ... TODO
...
'XGrid', 'off', ... [on | {off}]
'YGrid', 'off', ... [on | {off}]
'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
'XMinorGrid', 'off' , ... [on | {off}]
'YMinorGrid', 'off', ... [on | {off}]
'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
...
...'XTick', 0:0.1:100, ... ticks of x axis
...'YTick', 0:1:10, ... ticks of y axis
...'XTickLabel', {'-1','0','1'}, ...
...'YTickLabel', {'-1','0','1'}, ...
'XMinorTick', 'off', ... [on | {off}]
'YMinorTick', 'off', ... [on | {off}]
'TickDir', 'out', ... [{in} | out] inside or outside (for 2D)
'TickLength', [.01 .01], ... length of the ticks
...
'XColor', [.1 .1 .1], ... color of x axis
'YColor', [.1 .1 .1], ... color of y axis
'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
...'XLim', [0 1], ... limits for the x-axis
...'YLim', [0 1], ... limits for the y-axis
...
'FontName', strFontName, ... kind of fonts of labels
'FontSize', iFontSize, ... size of fonts of labels
'FontUnits', strFontUnit, ... units of the size of fonts
'FontWeight', strFontWeight, ... weight of fonts of labels
'FontAngle', strFontAngle, ... inclination of fonts of labels
...
'LineWidth', fLineWidth); % width of the line of the axes
% fonts properties
xlab = xlabel('z [ m ]');
set(xlab,'interpreter','Latex','FontSize',20)
ylab = ylabel('$u_{z} \times B_{y} [ V m^{-1} ]$');
set(ylab,'interpreter','Latex','FontSize',14)
subplot(2,2,4)
plot(z,j,'LineWidth',2)
axis tight
set(gca, ...
... 'Position', [1 1 20 10], ... TODO
... 'OuterPosition', [1 1 20 10], ... TODO
...
'XGrid', 'off', ... [on | {off}]
'YGrid', 'off', ... [on | {off}]
'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
'XMinorGrid', 'off' , ... [on | {off}]
'YMinorGrid', 'off', ... [on | {off}]
'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
...
...'XTick', 0:0.1:100, ... ticks of x axis
...'YTick', 0:1:10, ... ticks of y axis
...'XTickLabel', {'-1','0','1'}, ...
...'YTickLabel', {'-1','0','1'}, ...
'XMinorTick', 'off', ... [on | {off}]
'YMinorTick', 'off', ... [on | {off}]
'TickDir', 'out', ... [{in} | out] inside or outside (for 2D)
'TickLength', [.01 .01], ... length of the ticks
...
'XColor', [.1 .1 .1], ... color of x axis
'YColor', [.1 .1 .1], ... color of y axis
'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
...'XLim', [0 1], ... limits for the x-axis
...'YLim', [0 1], ... limits for the y-axis
...
'FontName', strFontName, ... kind of fonts of labels
'FontSize', iFontSize, ... size of fonts of labels
'FontUnits', strFontUnit, ... units of the size of fonts
'FontWeight', strFontWeight, ... weight of fonts of labels
'FontAngle', strFontAngle, ... inclination of fonts of labels
...
'LineWidth', fLineWidth); % width of the line of the axes
% fonts properties
xlab = xlabel('z [ m ]');
set(xlab,'interpreter','Latex','FontSize',20)
ylab = ylabel('$j_x[ A m^{-2} ]$');
set(ylab,'interpreter','Latex','FontSize',14)