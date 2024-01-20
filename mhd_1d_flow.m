clc; close all; clear all;
format shortE

%% Initial Parameters

% Global
% will need gravity
gamma = 5/3; % (Constant) Specific heat ratio for monoatomic gas
a = 340; % Sound speed [m/s]
R = 8.314; % Ideal gas constant [J/(mol*K)]
k_b = 8.617333262145e-5; % Boltzmann constant [eV/K]
ln_lambda = 10; % Coulomb logarithm (from paper)
alpha = 1e-3; % Loss of heat by current
muo = (4*pi)*10^(-7); % Permeability of a vaccuum [Henry/m]

% Engine
I = 3500; % Current [Amps]
mdot = 0.02; % Mass flow rate [g/s] 0.2
M = 1.1; % Mach number
T = 1; % Temperature [eV]
Tk = T/(k_b); % Temperature [Kelvin]
v_o = M*a; % Initial velocity [m/s]

% Geometry [m]
r_c = 0.02; % Cathode radius
r_o = 0.07; % Inlet anode radius
% Ly = pi*(r_c+r_o);
L = 0.112; % Length of duct
z = linspace(0, L, 1000).'; % Axial coordinate
A = (r_o+r_c)*L; % Area (r_a+r_o+r_c)*L; r_o+r_c

% Magnetic Field [Tesla]
Bo = 5e-2; % Initial magnetic field at Inlet
% B_self = Bo; % Self-induced magnetic field

% Electric Field
V = 10; % Supplied Voltage [Volts]
E = V/r_o; % [Volts/m]
sigma = ((1.5085e-2)*T^(3/2))/(ln_lambda);

%% Constant Parameters

rho = 1000/(v_o*A); % Density (constant) [g/m^3]
p_o = rho*R*Tk; % Initial Pressure [Pa] OR [J/m^3] at inlet
j_o = I/A; % Current density [Amps/m^2] sigma*(E-(v*B_self))

%% Solving Flow

Y0 = [v_o; p_o; Tk; j_o; Bo;]; % Inlet Boundary Conditions

[z_coordinates, solved_ode] = ode45(@(t, y) fun(Y0, ...
                                                sigma, ...
                                                rho, ...
                                                E, ...
                                                R, ...
                                                alpha, ...
                                                muo), ...
                                                z, ...
                                                Y0); % ODE Solver (Runge-Kutta)

% Normalised Parameters
U_star = solved_ode(:,1);
p_star = solved_ode(:,2);

% figure;
subplot(2,1,1)
plot(z_coordinates, U_star);
xlabel('Length of Thruster [m]');
ylabel('Velocity [m/s]');
title('Velocity over Thruster Length');

subplot(2,1,2)
plot(z_coordinates, p_star);
xlabel('Length of Thruster [m]');
ylabel('Pressure [Pa]');
title('Pressure Change over Thruster Length');

%% Functions

function dY_dz = fun(Y0, sigma, rho, E, R, alpha, muo)
    c_v = (3/2)*R; % Specific heat at specific volume (ideal monatomic)
    c_p = c_v+R; % Specific heat ratio at specific pressure (ideal monatomic)
    
    % (Inlet) z=0 Boundary Initial Conditions
    v = Y0(1);
    p = Y0(2);
    T = Y0(3);
    j = Y0(4);
    B_self = Y0(5);

    % F(Y) Right-Hand Side
    rhs_eq = zeros(size(Y0));
    rhs_eq(1) = (1/rho)*j*B_self;                          % (1) Momentum Equation
    rhs_eq(2) = (j*v*B_self)+(((j^2)/sigma)*(1*alpha));    % (2) Energy Equation 
    rhs_eq(3) = 0;                                         % (3) Equation of State rho*R*T
    rhs_eq(4) = 0;                                         % (4) Ohm's Law (Current Distribution) sigma*(E-(v*B_self));
    rhs_eq(5) = -muo*j;                                    % (5) Maxwell's equation for self induced mag field

%     disp(rhs_eq)

    % Coninuity equation equation taken out as no density diff term involved.
    % Mass matrix missing r/r_o - uniform cross secton so will come back

    % Mass Matrix
    M = zeros(length(Y0), length(Y0)); % 5x5
    M(1,1) = v;             M(1,2) = (1/rho);                   % (1) Momentum Equation
    M(2,1) = rho*(v^2);     M(2,4) = rho*v*c_p;                 % (2) Energy Equation
    M(3,2) = 1;             M(3,4) = -(rho*R);                  % (3) Equation of State
    M(4,1) = B_self*sigma;  M(4,3) = sigma*v;   M(4,5) = 1;     % (4) Ohm's Law (Current Distribution)
    M(5,3) = 1;                                                 % (5) B_self

%     disp(M)

    dY_dz = M \ rhs_eq;
    disp(dY_dz)
end











