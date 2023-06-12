% v0.1 RadCom simulator
% @Author: Marco Manzoni (marco.manzoni@polimi.it)

clear all;
close all;
clc;

%% Antenna system parameter
f0                  = 4e9;              % Carrier frequency [Hz]
c                   = physconst("lightspeed");
K                   = physconst('boltzmann'); % Boltzmann constant
lambda              = c/f0;

Nx_tx               = 1;                % Number of antennae in the horizontal plane
Nz_tx               = 1;                % Number of antennae in the vertical plane
dx                  = lambda/4;         % Antenna spacing in the horizontal plane
dz                  = dx;               % Antenna spacing in the vertical plane
installation_height = 20;               % The height of the antenna over the ground plane
delta_psi           = 100/180*pi;       % azimuth beamwidth [rad]
delta_teta          = 100/180*pi;       % elevation beamwidth [rad]
teta_point          = 0/180*pi;         % antenna pointing in elevation [rad] 
psi_point           = 0/180*pi;         % antenna pointing in azimuth [rad]

Ptx_db              = -30;              % Total transmission power [dBw] (dBw = dBm - 30)
eff_ant             = 0.7;              % Antenna efficiency
T0_ant              = 290;              % Antenna temperature [K]
T_scene             = 290;              % Scene temperature [K]
F_dB                = 10;                % Receiver noise figure [dB]
T0                  = 290;              % Temperature at which the noise figure is referring to [K]

%% Parameters for the 5G sub-6 system

B                 = 30e6;               % Bandwidth [Hz]
u                 = 0;                  % Numerology

delta_f           = 15e3 * 2^u;         % Sub-carrier spacing (defined in this way by the standard)
Ts_no_cp          = 1/delta_f;          % Length of the OFDM symbol [s]

cpPercentage      = 7;                  % Percentage of the symbol dedicated to cyclic prefix
Ts_tot            = Ts_no_cp*(1+cpPercentage/100); % Total symbol time, including cyclic prefix

subCarrierPerRB   = 12;                 % Number of sub-carrier for each resource block;
numberOfRB        = ceil(B/(subCarrierPerRB*delta_f)); % Total number of RB in frequency;
totalSubcarriers  = ceil(B/delta_f);    % Total number of sub-carriers
frameDuration     = 10e-3;              % Duration of the frame [s];
Nsymbols          = floor(frameDuration/Ts_tot); 

subframePerFrame  = 10;                 % Number of sub-frame in the frame
subFrameDuration  = frameDuration/subframePerFrame; % Duration of each sub-frame [s]

slotsPerSubFrame  = 2^u;                % Number of slots for each sub frame;
slotDuration      = subFrameDuration/slotsPerSubFrame; % Duration of the slot [s];

symbolsPerSlot    = 14;                 % Number of symbols for each time slot
symbolDuration    = slotDuration/symbolsPerSlot; % Duration of each symbol [s];

PRF               = 1/symbolDuration;   % Pulse repetition frequency [Hz]
M                 = 8;                  % QPSK signal constellation size

% Generation of the data
frame = randsrc(totalSubcarriers, Nsymbols, 0:M-1);
figure; imagesc(frame); colorbar; xlabel("OFDM symbol number"); ylabel("Sub-carrier number");
title("OFDM Frame"); axis xy tight

% Associate to each number of QPSK symbol
%qpsk_modulated_data = qammod(data_source, M, UnitAveragePower=true);
qpsk_modulated_data = pskmod(frame, M);

% Simulate continuous time signal
fs                  = 2*B;     % Sampling frequency in fast time [Hz];
dt                  = 1/fs;   % sampling in fast time [s]
t                   = -Ts_tot:dt:Ts_tot;% Fast time axis
r                   = t*c/2; % Range axis
rho_rg              = c/(2*B); % Slant range resolution c/2B

% Define the frequency domain axis
NFFT                = 2048;
df                  = fs/NFFT;
f_ax                = (-NFFT/2:NFFT/2-1)*df;

% Along the first dimension I have the signal transmitted for a symbol in a
% sub-carrier, in the second simension I sweep the symbols (same
% sub-carrier) and along the third dimension I sweep all the sub-carrier.
% If I sum along the third dimension I get for every colum the OFDM symbol.
signal_tx = zeros(length(t), Nsymbols, totalSubcarriers);

for ii = 1:totalSubcarriers
    signal_tx(:,:,ii) = qpsk_modulated_data(ii,:).*exp(1j.*2*pi*(ii-1)*delta_f*t').*rectpuls(t'/Ts_tot);
end

% Visualize the spectrum of all the sub-carriers;
figure; plot(f_ax, db(fftshift(fft(squeeze(signal_tx(:,1,:)),NFFT)))); grid on; hold on;
xlabel("Frequency [Hz]"); title("Subcarriers of first OFDM symbol");

% Visualize the sum
plot(f_ax, db(fftshift(fft(sum(signal_tx(:,1,:),3),NFFT))), 'r', LineWidth=2); 

% Sum all the sub-carriers and bring it to passband for transmission
% The recived signal will be the delayed and scaled version of this signal
signal_tx = sum(signal_tx, 3).*exp(1j*2*pi*f0*t(:));

%% Trajectory of the platform

Sx = (0:Nsymbols-1)*lambda/4;
Sy = zeros(size(Sx));
Sz = zeros(size(Sx));

%% Target parameters
RCS                 = 1; % Radar cross section [m^2]

% Each line is a target, the colums are the x,y,z position. The center of
% the reference system is the center of the base station
p_t                 = [0, 100, -installation_height]; 

%% Derive some parameters from numbers above

Nt                  = size(p_t, 1);             % Number of targets

beam_sector         = delta_teta*delta_psi;     % [rad^2]
G                   = eff_ant*4*pi/beam_sector; % Antenna gain

Ptx = 10^(Ptx_db/10);                           %[W]

% Antenna positions in reception. The transmission antenna is located at
% the origin and transmits all the power
x_s = (0:Nx_tx-1)'*dx; x_s = x_s - mean(x_s); x_s = kron(ones(Nz_tx,1), x_s);
z_s = (0:Nz_tx-1)'*dz; z_s = z_s - mean(z_s); z_s = kron(z_s, ones(Nx_tx,1));
y_s = zeros(size(x_s));

N_ant       = length(x_s);                      % Total number of antennas
Ae          = (lambda^2)/4/pi*G;                % Equivalent area for each antenna

F           = 10^(F_dB/10);                     % noise figure in linear unit (noise factor)
Trx         = (F-1)*T0;                         % Receiver temperature [Kelvin]
T_ant       = eff_ant*T_scene + (1-eff_ant)*T0_ant; % Antenna temperature [Kelvin]
T_sys       = T_ant + Trx;                      % System temperature [Kelvin]
N0          = K*T_sys;                          % Noise Power Spectral Density [W/Hz]

%% Acquisition

Draw    = zeros(length(t), N_ant, Nsymbols);
Draw_n  = zeros(length(t), N_ant, Nsymbols);
Drc     = zeros(length(t), N_ant, Nsymbols);
Drc_n   = zeros(length(t), N_ant, Nsymbols);

for kk = 1:Nsymbols
    for ii = 1:N_ant
        for n = 1:Nt

            % In transmission
            delta_x = Sx(kk)-p_t(n,1);
            delta_y = Sy(kk)-p_t(n,2);
            delta_z = Sz(kk)-p_t(n,3);

            R_tx = sqrt(delta_x^2+delta_y^2+delta_z^2);

            % Azimuth angle from the radar to the target
            psi = asin(delta_x/R_tx);

            % off nadir angle from the radar to the target
            teta = asin(delta_y/(R_tx*cos(psi)));

            % Antenna pattern
            f = sinc((psi-psi_point)/delta_psi)^2*sinc((teta-teta_point+pi/2)/delta_teta)^2;

            % Power density at the target [W/m^2]
            Si = Ptx./(4*pi*R_tx.^2)*G.*f;

            % At the receiver
            delta_x = (Sx(kk)+x_s(ii))-p_t(n,1);
            delta_y = (Sy(kk)+y_s(ii))-p_t(n,2);
            delta_z = (Sz(kk)+z_s(ii))-p_t(n,3);

            R_rx = sqrt(delta_x^2+delta_y^2+delta_z^2);

            % Scattered power density at the receiver from each target  [W/m^2]
            Ss = Si./(4*pi*R_rx.^2).*RCS(n);

            % Power at the receiver from each target [W]
            Prx = Ss*Ae.*f; % f should be computed again for each receiver... approximation

            tau             = (R_tx+R_rx)/c; % Delay from distance (to change in case of bistatic)
            delayFilter     = fs*sinc(fs*(t-tau)); % Filter that I use for delay
            s               = sqrt(Prx)*conv(signal_tx(:,kk), delayFilter, 'same')*dt; % Delay the TX signal
            s               = s.*exp(-1j*2*pi*f*t(:)); % Demodulate it
            Draw(:,ii,kk)   = Draw(:,ii,kk) + s; % Sum for every target
        end

        Pn                  = N0*fs; % noise power, notice that fs is the simulation bandwidth
        n                   = sqrt(Pn/2)*(randn(length(t),1) + 1j.*randn(length(t),1)); % I generate noise
        Draw_n(:,ii, kk)    = Draw(:,ii, kk) + n; % and add it to the raw data

        % Range compression for the ii-th antenna
        Drc(:, ii, kk)      = conv(Draw(:,ii,kk), flip(conj(signal_tx(:,kk))), "same")*dt;
        Drc_n(:, ii, kk)    = conv(Draw_n(:,ii,kk), flip(conj(signal_tx(:,kk))), "same")*dt;
    end
end

figure; plot(r, db(Drc(:,:,1))); grid on; xlabel("Range [m]"); ylabel("Amplitude");
title("Range compressed data for each antenna without noise");

figure; plot(r, db(Drc_n(:,:,1))); grid on; xlabel("Range [m]"); ylabel("Amplitude");
title("Range compressed data for each antenna with noise");

sqrt(Prx)*(Ts_tot*totalSubcarriers) % Theoretical peak value that I should get
mean(max(abs(Drc(:,:,1)))) % Experimental peak value that I get

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


































