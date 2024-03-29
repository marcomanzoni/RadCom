% v0.2 RadCom simulator
% @Author: Marco Manzoni (marco.manzoni@polimi.it)

clear all;
close all;
clc;

load('./lajolla.mat');
lajolla = flip(lajolla);

addpath("C:\Users\manzoni\OneDrive - Politecnico di Milano\Progetti\JC&S\codice");

%% Antenna system parameter
f0                  = 5.9e9;              % Carrier frequency [Hz]
c                   = physconst("lightspeed");
K                   = physconst('boltzmann'); % Boltzmann constant
lambda              = c/f0;

Nx_tx               = 1;                % Number of antennae in the horizontal plane
Nz_tx               = 1;                % Number of antennae in the vertical plane
dx                  = lambda/4;         % Antenna spacing in the horizontal plane
dz                  = dx;               % Antenna spacing in the vertical plane
installation_height = 50;               % The height of the antenna over the ground plane
delta_psi           = 100/180*pi;       % azimuth beamwidth [rad]
delta_teta          = 70/180*pi;       % elevation beamwidth [rad]
teta_point          = 45/180*pi;         % antenna pointing in elevation [rad]
psi_point           = 0/180*pi;         % antenna pointing in azimuth [rad]

Ptx_db              = -24;              % Total transmission power [dBw] (dBw = dBm - 30)
eff_ant             = 0.7;              % Antenna efficiency
T0_ant              = 290;              % Antenna temperature [K]
T_scene             = 290;              % Scene temperature [K]
F_dB                = 7;               % Receiver noise figure [dB]
T0                  = 290;              % Temperature at which the noise figure is referring to [K]

%% Parameters for the 5G sub-6 system

B                 = 40e6;               % Bandwidth [Hz]
u                 = 3;                  % Numerology

delta_f           = 15e3 * 2^u;         % Sub-carrier spacing (defined in this way by the standard)
Ts_no_cp          = 1/delta_f;          % Length of the OFDM symbol [s]

cpPercentage      = 0;                  % Percentage of the symbol dedicated to cyclic prefix
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
M                 = 4;                  % QPSK signal constellation size

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
rho_az = rho_rg;

% Define the frequency domain axis to check the bandwidth occupancy
NFFT                = 2^nextpow2(length(t));
df                  = fs/NFFT;
f_ax                = (-NFFT/2:NFFT/2-1)*df;

% Along the first dimension I have the signal transmitted for a symbol in a
% sub-carrier, in the second simension I sweep the symbols (same
% sub-carrier) and along the third dimension I sweep all the sub-carrier.
% If I sum along the third dimension I get for every colum the OFDM symbol.
subCarrierIndex = (0:totalSubcarriers-1); subCarrierIndex = floor(subCarrierIndex - mean(subCarrierIndex));
modulationMatrix = exp(1j.*2*pi*t'*subCarrierIndex*delta_f);
signal_tx_bb = modulationMatrix*qpsk_modulated_data.*rectpuls(t'/Ts_tot);

% Visualize the spectrum of all the sub-carriers;
figure; plot(f_ax, db(fftshift(fft(squeeze(signal_tx_bb(:,1:5)),NFFT)))); grid on; hold on;
xlabel("Frequency [Hz]"); title("First 5 OFDM symbols");

y = autocorrelation(signal_tx_bb(:,1), signal_tx_bb(:,1), dt);
figure; plot(r,db(y)); grid on; xlim([-100,100]); ylabel("Amplitude [dB]"); xlabel("Range [m]");
title("Autocorrelation of the first transmitted symbol");

% The recived signal will be the delayed and scaled version of this signal.
% First we have to scale it to have unit power
Eg = sum(abs(signal_tx_bb).^2)*dt;
Pw = Eg/Ts_tot;
signal_tx = signal_tx_bb./sqrt(Pw);%.*exp(1j*2*pi*f0*t(:));

%% Target parameters
RCS                 = 1; % Radar cross section [m^2]

% Each line is a target, the colums are the x,y,z position. The center of
% the reference system is the center of the base station. If you place the
% target below the ground it is intended under the snow

snow_depth = 2;

%p_t                 = [0, 150, -installation_height]; %On the ground
p_t                 = [0, 150, -installation_height-snow_depth]; %three meters under the snow

R_t = sqrt(p_t(1).^2 + p_t(2).^2 + p_t(3).^2);

%% Trajectory of the platform
As = lambda/2/rho_az*R_t;
vp = 2; % [m/s]

Sx = -As/2 : vp/PRF :As/2;% linspace(-As/2, As/2, Nsymbols);
Sy = 0*Sx;
Sz = 0*Sx;

Ntau = length(Sx)

%% Derive some parameters from numbers above

Nt                  = size(p_t, 1);             % Number of targets
beam_sector         = delta_teta*delta_psi;     % [rad^2]
G                   = eff_ant*4*pi/beam_sector; % Antenna gain

Ptx                 = 10^(Ptx_db/10);               %[W]
EIRP_w              = Ptx*G;                        %[W]
EIRP_dbm            = 10*log10(EIRP_w*1000)         % [dBm]

% Antenna positions in reception. The transmission antenna is located at
% the origin and transmits all the power
x_s = (0:Nx_tx-1)'*dx; x_s = x_s - mean(x_s); x_s = kron(ones(Nz_tx,1), x_s);
z_s = (0:Nz_tx-1)'*dz; z_s = z_s - mean(z_s); z_s = kron(z_s, ones(Nx_tx,1));
y_s = zeros(size(x_s));

N_ant       = length(x_s);                          % Total number of antennas
Ae          = (lambda^2)/4/pi*G;                    % Equivalent area for each antenna

F           = 10^(F_dB/10);                         % noise figure in linear unit (noise factor)
Trx         = (F-1)*T0;                             % Receiver temperature [Kelvin]
T_ant       = eff_ant*T_scene + (1-eff_ant)*T0_ant; % Antenna temperature [Kelvin]
T_sys       = T_ant + Trx;                          % System temperature [Kelvin]
N0          = K*T_sys;                              % Noise Power Spectral Density [W/Hz]

% Calculat the one-way attenuation constant in the snow
W = 0.005;
alpha       = snowPowerAttenuationLuca(f0,W);

%% Acquisition

Draw    = zeros(length(t), N_ant, Ntau, 'single');
Draw_n  = zeros(length(t), N_ant, Ntau, 'single');
Drc     = zeros(length(t), N_ant, Ntau, 'single');
Drc_n   = zeros(length(t), N_ant, Ntau, 'single');

idx_symbol = 1;

for kk = 1:Ntau % slow time samples
    for ii = 1:N_ant % Antennae
        for n = 1:Nt % Targets

            fprintf("%d / %d \n", kk, Ntau);

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

            % Calculate the approximate path length under the snow;
            l = snow_depth;%pathLengthUnderSnow(Sx(kk), Sy(kk), Sz(kk), p_t(n,1), p_t(n,2), p_t(n,3), installation_height);

            L = exp(-2*alpha*l);

            % Power density at the target [W/m^2]
            Si = Ptx./(4*pi*R_tx.^2)*G.*f*L;

            % At the receiver
            delta_x = (Sx(kk)+x_s(ii))-p_t(n,1);
            delta_y = (Sy(kk)+y_s(ii))-p_t(n,2);
            delta_z = (Sz(kk)+z_s(ii))-p_t(n,3);

            R_rx = sqrt(delta_x^2+delta_y^2+delta_z^2);

            % Scattered power density at the receiver from each target  [W/m^2]
            Ss = Si./(4*pi*R_rx.^2).*RCS(n)*L;

            % Power at the receiver from each target [W]
            Prx = Ss*Ae.*f; % f should be computed again for each receiver... approximation

            tau             = (R_tx+R_rx)/c; % Delay from distance (to change in case of bistatic)

            s = sqrt(Prx)*interp1(t', signal_tx(:,idx_symbol), t'-tau, "cubic", 0).*exp(-1j*2*pi*f0*tau);
            %y = autocorrelation(s, signal_tx(:,kk), dt);
            %figure; plot(r, abs(y));

            %delayFilter     = fs*sinc(fs*(t-tau)); % Filter that I use for delay
            %s               = sqrt(Prx)*conv(signal_tx(:,kk), delayFilter, 'same')*dt; % Delay the TX signal
            %s               = s.*exp(-1j*2*pi*f*t(:)); % Demodulate it
            
            Draw(:,ii,kk)   = Draw(:,ii,kk) + s; % Sum for every target
        end

        Pn                  = N0*fs; % noise power, notice that fs is the simulation bandwidth
        n                   = sqrt(Pn/2)*(randn(length(t),1) + 1j.*randn(length(t),1)); % I generate noise
        Draw_n(:,ii, kk)    = Draw(:,ii, kk) + n; % and add it to the raw data

        % Range compression for the ii-th antenna
        Drc(:, ii, kk)      = conv(Draw(:,ii,kk), flip(conj(signal_tx(:,idx_symbol))), "same")*dt;
        Drc_n(:, ii, kk)    = conv(Draw_n(:,ii,kk), flip(conj(signal_tx(:,idx_symbol))), "same")*dt;

        idx_symbol = idx_symbol + 1;
        if idx_symbol > Nsymbols
            idx_symbol = 1;
        end
    end
end

figure; plot(r, db(Drc(:,:,1))); grid on; xlabel("Range [m]"); ylabel("Amplitude");
title("Range compressed data for each antenna without noise");

figure; plot(r, db(Drc_n(:,:,1))); grid on; xlabel("Range [m]"); ylabel("Amplitude");
title("Range compressed data for each antenna with noise");

sqrt(Prx)*(Ts_tot) % Theoretical peak value that I should get
mean(max(abs(Drc(:,:,1)))) % Experimental peak value that I get

%% Image focusing

Nres = 3;

x_foc = min(p_t(:,1))-Nres*rho_rg : 0.05*rho_rg : max(p_t(:,1))+Nres*rho_rg;
y_foc = min(p_t(:,2))-Nres*rho_rg : 0.05*rho_rg : max(p_t(:,2))+Nres*rho_rg;

[Y,X] = ndgrid(y_foc,x_foc);
Z = p_t(3)*ones(size(X)); % I decide to focus on ground, but since I have 3D resolution I could focus other planes

I = zeros(size(X), 'single');
I_n = zeros(size(I));

for kk = 1:Ntau % symbols = slow time samples

    fprintf("Focusing slow-time %d/%d: ", kk, Ntau);

    delta_x = Sx(kk)-X;
    delta_y = Sy(kk)-Y;
    delta_z = Sz(kk)-Z;

    R_tx = sqrt(delta_x.^2+delta_y.^2+delta_z.^2);

    for ii = 1:N_ant
        delta_x = (Sx(kk)+x_s(ii))-X;
        delta_y = (Sy(kk)+y_s(ii))-Y;
        delta_z = (Sz(kk)+z_s(ii))-Z;

        R_rx = sqrt(delta_x.^2+delta_y.^2+delta_z.^2);

        tau   = (R_tx + R_rx)/c; % Delay from distance (to change in case of bistatic)

        I = I + interp1(t, Drc(:, ii, kk) , tau,"linear",NaN).*exp(+1j.*2*pi*f0*tau);
        I_n = I_n + interp1(t, Drc_n(:,ii,kk), tau,"linear",NaN).*exp(+1j.*2*pi*f0*tau);
    end
    fprintf("Done. \n");
end

% figure; imagesc(x_foc, y_foc, abs(I(:,:,1))); axis xy tight; title("MIMO Antenna without noise");
% xlabel("x [m]"); ylabel("y [m]"); colorbar; colormap(lajolla);
% 
% figure; imagesc(x_foc, y_foc, abs(I_n(:,:,1))); axis xy tight; title("MIMO Antenna with noise");
% xlabel("x [m]"); ylabel("y [m]"); colorbar; colormap(lajolla);


figure; imagesc(x_foc, y_foc, abs(I)); axis xy tight; title({sprintf("Depth under snow = %.2f m", snow_depth), sprintf("UAV altitude = %.2f m", installation_height)});
xlabel("x [m]"); ylabel("y [m]"); colormap(lajolla);

figure; imagesc(x_foc, y_foc, abs(I_n)); axis xy tight; title({sprintf("Depth under snow = %.2f m", snow_depth), sprintf("UAV altitude = %.2f m", installation_height)});
xlabel("x [m]"); ylabel("y [m]"); colormap(lajolla); axis equal tight
exportgraphics(gca, fullfile("C:\Users\manzoni\OneDrive - Politecnico di Milano\Conferenze\[2024] EuSAR\Immagini", sprintf("SAR_noise_%.2f_%d.pdf", snow_depth, installation_height)));

































