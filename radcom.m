% v0.1 RadCom simulator
% @Author: Marco Manzoni (marco.manzoni@polimi.it)

clear all;
close all;
clc;

%% Antenna system parameter
f0                  = 4e9; % Carrier frequency [Hz]
c                   = physconst("lightspeed");
K                   = physconst('boltzmann'); % Boltzmann constant
lambda              = c/f0;

Nx_tx               = 8;                % Number of antennae in the horizontal plane
Nz_tx               = 4;                % Number of antennae in the vertical plane
dx                  = lambda/4;         % Antenna spacing in the horizontal plane
dz                  = dx;               % Antenna spacing in the vertical plane
installation_height = 20;               % The height of the antenna over the ground plane
delta_psi           = 100/180*pi;       % azimuth beamwidth [rad]
delta_teta          = 100/180*pi;       % elevation beamwidth [rad]
teta_point          = 0/180*pi;        % antenna pointing in elevation [rad] 
psi_point           = 0/180*pi;         % antenna pointing in azimuth [rad]

Ptx_db              = -30;                % Total transmission power [dBw]
eff_ant             = 0.7;              % Antenna efficiency
T0_ant              = 290;              % Antenna temperature [K]
T_scene             = 290;              % Scene temperature [K]
F_dB                = 3;                % Receiver noise figure [dB]
T0                  = 290;              % Temperature at which the noise figure is referring to [K]

%% OFDM parameters
Ts                = 12e-6;                             % Length of the OFDM symbol [s]
M                 = 2;                  % QPSK signal constellation size
no_of_data_points = 128;                % Number of symbols going out from the source
block_size        = 128;                % size of each ofdm block = number of subcarriers
no_of_ifft_points = block_size;         % Points for the FFT/IFFT
no_of_fft_points  = block_size;

%% Target parameters
RCS                 = 1; % Radar cross section [m^2]

% Each line is a target, the colums are the x,y,z position. The center of
% the reference system is the center of the base station
p_t                 = [0, 100, -installation_height]; 

%% Derive some parameters from numbers above

Nt                  = size(p_t, 1); % Number of targets

beam_sector         = delta_teta*delta_psi; % [rad^2]
G                   = eff_ant*4*pi/beam_sector; % Antenna gain

Ptx = 10^(Ptx_db/10);

% Antenna positions
x_s = (0:Nx_tx-1)'*dx; x_s = x_s - mean(x_s); x_s = kron(ones(Nz_tx,1), x_s);
z_s = (0:Nz_tx-1)'*dz; z_s = z_s - mean(z_s); z_s = kron(z_s, ones(Nx_tx,1));
y_s = zeros(size(x_s));

N_ant       = length(x_s);          % Total number of antennas
P_tx_ant    = Ptx/N_ant;            % Power for each antenna
Ae          = (lambda^2)/4/pi*G;    % Equivalent area for each antenna

F           = 10^(F_dB/10); % noise figure in linear unit (noise factor)
Trx         = (F-1)*T0; % Receiver temperature [Kelvin]
T_ant       = eff_ant*T_scene + (1-eff_ant)*T0_ant; % Antenna temperature [Kelvin]
T_sys       = T_ant + Trx; % System temperature [Kelvin]
N0          = K*T_sys; % Noise Power Spectral Density [W/Hz]

figure;
plot3(x_s, y_s, z_s, 'rd'); grid on; hold on;
plot3(p_t(:,1), p_t(:,2), p_t(:,3), 'b*');
xlabel("x [m]"); ylabel("y [m]"); zlabel("z [m]");

%% Generation of the symbols to be transmitted

% Generate "no_of_data_points" random numbers from 0 to M-1
data_source = randsrc(1, no_of_data_points, 0:M-1);

% Associate to each number of QPSK symbol
%qpsk_modulated_data = qammod(data_source, M, UnitAveragePower=true);
qpsk_modulated_data = pskmod(data_source, M);

delta_f             = 1/Ts;     % Spacing between subcarriers
ovs                 = 512;      % Over Sampling Factor w.r.t. Ts to simulate continuous time
dt                  = Ts/ovs;   % sampling in fast time [s]
fs                  = 1/dt;     % Sampling frequency in fast time [Hz];
t                   = -Ts:dt:Ts;% Fast time axis
r                   = t*c/2; % Range axis
rho_rg              = c/(2*delta_f*block_size); % Slant range resolution c/2B

% In each column there is the signal transmitted in a subcarrier, the final
% transmitted signal is the sum of all the signals in all the sub-carriers.
signal_tx = zeros(length(t), block_size);
for ii = 1:block_size
    signal_tx(:,ii) = qpsk_modulated_data(ii).*exp(1j.*2*pi*(ii-1)*delta_f*t).*rectpuls(t/Ts);
end
% Sum all the sub-carriers and bring it to passband for transmission
signal_tx = sum(signal_tx, 2).*exp(1j*2*pi*f0*t(:));

% The recived signal will be the delayed and scaled version of this signal

%% Acquisition

Draw    = zeros(length(t), N_ant);
Draw_n  = zeros(length(t), N_ant);
Drc     = zeros(length(t), N_ant);
Drc_n   = zeros(length(t), N_ant);

for ii = 1:N_ant
     for n = 1:Nt
        delta_x = x_s(ii)-p_t(n,1);
        delta_y = y_s(ii)-p_t(n,2);
        delta_z = z_s(ii)-p_t(n,3);

        R = sqrt(delta_x^2+delta_y^2+delta_z^2);

        % Azimuth angle from the radar to the target
        psi = asin(delta_x/R);

        % off nadir angle from the radar to the target
        teta = asin(delta_y/(R*cos(psi)));

        % Antenna pattern
        f = sinc((psi-psi_point)/delta_psi)^2*sinc((teta-teta_point+pi/2)/delta_teta)^2;

        % Power density at the target [W/m^2]
        Si = P_tx_ant./(4*pi*R.^2)*G.*f;

        % Scattered power density at the receiver from each target  [W/m^2]
        Ss = Si./(4*pi*R.^2).*RCS(n);

        % Power at the receiver from each target [W]
        Prx = Ss*Ae.*f;

        tau         = 2*R/c; % Delay from distance (to change in case of bistatic)
        delayFilter = fs*sinc(fs*(t-tau)); % Filter that I use for delay
        s           = sqrt(Prx)*conv(signal_tx, delayFilter, 'same')*dt; % Delay the TX signal
        s           = s.*exp(-1j*2*pi*f*t(:)); % Demodulate it
        Draw(:,ii)  = Draw(:,ii) + s; % Sum for every target 
     end

     Pn             = N0*fs; % noise power, notice that fs is the simulation bandwidth
     n              = sqrt(Pn/2)*(randn(length(t),1) + 1j.*randn(length(t),1)); % I generate noise
     Draw_n(:,ii)   = Draw(:,ii) + n; % and add it to the raw data

     % Range compression for the ii-th antenna
     Drc(:,ii)      = conv(Draw(:,ii), flip(conj(signal_tx)), "same")*dt;
     Drc_n(:,ii)    = conv(Draw_n(:,ii), flip(conj(signal_tx)), "same")*dt;
end

figure; plot(r, db(Drc)); grid on; xlabel("Range [m]"); ylabel("Amplitude");
title("Range compressed data for each antenna without noise");

figure; plot(r, db(Drc_n)); grid on; xlabel("Range [m]"); ylabel("Amplitude");
title("Range compressed data for each antenna with noise");

sqrt(Prx)*(Ts*block_size) % Theoretical peak value that I should get
mean(max(abs(Drc))) % Experimental peak value that I should get

%% Image focusing

x_foc = min(p_t(:,1))-30 : 0.1 : max(p_t(:,1))+30;
y_foc = min(p_t(:,2))-30 : 0.1 : max(p_t(:,2))+30;

[Y,X] = ndgrid(y_foc,x_foc);
Z = zeros(size(X)); % I decide to focus on ground, but since I have 3D resolution I could focus other planes

I = zeros(size(X));
I_n = zeros(size(I));

for ii = 1:N_ant
    delta_x = x_s(ii)-X;
    delta_y = y_s(ii)-Y;
    delta_z = z_s(ii)-Z;

    R     = sqrt(delta_x.^2+delta_y.^2+delta_z.^2);
    tau   = 2*R/c; % Delay from distance (to change in case of bistatic)

    I = I + interp1(t, Drc(:,ii), tau,"linear",NaN).*exp(+1j.*2*pi*f0*tau);
    I_n = I_n + interp1(t, Drc_n(:,ii), tau,"linear",NaN).*exp(+1j.*2*pi*f0*tau);
end

figure; imagesc(x_foc, y_foc, abs(I)); axis xy tight; title("Focused image without noise");
xlabel("x [m]"); ylabel("y [m]");

figure; imagesc(x_foc, y_foc, abs(I_n)); axis xy tight; title("Focused image with noise");
xlabel("x [m]"); ylabel("y [m]");


































