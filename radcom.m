% v0.1 RadCom simulator
% @Author: Marco Manzoni (marco.manzoni@polimi.it)

clear all;
close all;
clc;

%% Antenna system parameter
f0                  = 4e9; % Carrier frequency [Hz]
c                   = physconst("lightspeed");
lambda              = c/f0;

Nx_tx               = 8; % Number of antennae in the horizontal plane
Ny_tx               = 4; % Number of antennae in the vertical plane
dx                  = lambda/4; % Antenna spacing in the horizontal plane
dy                  = dx; % Antenna spacing in the vertical plane
installation_height = 20; % The height of the antenna over the ground plane
delta_psi           = 100/180*pi; % azimuth beamwidth [rad]
delta_teta          = 100/180*pi; % elevation beamwidth [rad]
teta_point          = 0/180*pi; % antenna pointing in elevation [rad] 
psi_point           = 0/180*pi;

Ptx_db              = 0; % Total transmission power [dBm]
eff_ant             = 0.7; % Antenna efficiency
T0_ant              = 290; % Antenna temperature [K]
F                   = 3; % Receiver noise figure [dB]
T0                  = 290; % Temperature at which the noise figure is referring to [K]

%% Target parameters
RCS                 = 10; % Radar cross section [m^2]

% Each line is a target, the colums are the x,y,z position. The center of
% the reference system is the center of the base station
p_t                 = [0, -installation_height, 20]; 

%% Derive some parameters from numbers above
Nt                  = size(p_t, 1);

beam_sector         = delta_teta*delta_psi; % [rad^2]
G                   = eff_ant*4*pi/beam_sector; % Antenna gain

Ptx = 10^(Ptx_db/10);

% Antenna positions
x_s = (0:Nx_tx-1)'*dx; x_s = x_s - mean(x_s); x_s = kron(ones(Ny_tx,1), x_s);
y_s = (0:Ny_tx-1)'*dy; y_s = y_s - mean(y_s); y_s = kron(y_s, ones(Nx_tx,1));
N_ant = length(x_s);

%% Acquisition

for n = 1:N
     for ii = 1:N_ant
        delta_x = x_s(ii)-p_t(n,1);
        delta_y = y_s(ii)-p_t(n,2);
        delta_z = z_s(ii)-p_t(n,3);

        R = sqrt(delta_x^2+delta_y^2+delta_z^2);


     end
end
