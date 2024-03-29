clear all;
%close all;
clc;

%% Antenna system parameter
f0                  = 5.9e9;              % Carrier frequency [Hz]
c                   = physconst("lightspeed");
K                   = physconst('boltzmann'); % Boltzmann constant
lambda              = c/f0;

%% Parameters for the 5G sub-6 system

u_vect = 0:4;

for kk = 1:length(u_vect)

    B                 = 40e6;               % Bandwidth [Hz]
    u                 = u_vect(kk);                  % Numerology
    rho_rg            = c/2/B;

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
    M                 = 4;                  % QAM signal constellation size. !!!!!!!!!CHANGE THIS!!!!!!

    % Generation of the data
    frame = randsrc(totalSubcarriers, Nsymbols, 0:M-1);

    % Associate to each number of QPSK symbol
    qpsk_modulated_data = qammod(frame, M, UnitAveragePower=true);
    %qpsk_modulated_data = pskmod(frame, M);

    % Simulate continuous time signal
    fs                  = 30*B;     % Sampling frequency in fast time [Hz];
    dt                  = 1/fs;   % sampling in fast time [s]
    t                   = -Ts_tot:dt:Ts_tot;% Fast time axis
    r                   = t*c/2; % Range axis
    dr = mean(diff(r));

    % Along the first dimension I have the signal transmitted for a symbol in a
    % sub-carrier, in the second simension I sweep the symbols (same
    % sub-carrier) and along the third dimension I sweep all the sub-carrier.
    % If I sum along the third dimension I get for every colum the OFDM symbol.
    subCarrierIndex = (0:totalSubcarriers-1); subCarrierIndex = floor(subCarrierIndex - mean(subCarrierIndex));
    modulationMatrix = exp(1j.*2*pi*t'*subCarrierIndex*delta_f);
    signal_tx_bb = modulationMatrix*qpsk_modulated_data.*rectpuls(t'/Ts_tot);

    %% Pulse compression

    clc;

    N = 100;

    % I take only one OFDM symbol
    s  = signal_tx_bb(:,1:N);

    % Normalize its power
    Es = sum(abs(s).^2)*dt;
    s  = s./sqrt(Es/Ts_tot);

    % Number of freq points
    Nf = 2^nextpow2(length(t));

    % Frequency analysis
    [S,f_ax] = dft(s,t,Nf); % Take the DFT
    S        = S*dt; % Adjust its amplitude
    SS       = S.*conj(S); % Perform autocorrelation in freq domain
    Drc_f    = idft(SS,f_ax,t); % Back in time domain
    df       = f_ax(2)-f_ax(1);
    Drc_f    = Drc_f*df*Nf; % Adjust once again the amplitudes

    % Pulse compression with convolution in time domain
    Drc_t = zeros(size(s));
    for ii = 1:size(s,2)
        Drc_t(:,ii) = conv(s(:,ii),flipud(conj(s(:,ii))),'same')*dt;
    end
    %Drc_t = mean(abs(Drc_t),2);

    [~, indexRes] = min(abs(r-rho_rg));
    [~, indexZero] = min(abs(r));
    delta = indexRes-indexZero;

    % Main lobe I say 3 resolution cells
    indMainLobe = indexZero-2*delta : indexZero+2*delta;

    %
    indTenRes = indexZero-10*delta : indexZero+10*delta;

    % Everything else, apart the main lobe
    indOutside = [1:indMainLobe(1)-1, indMainLobe(end)+1:length(Drc_t)];

    PSLR_time(kk)       = 20*log10(mean(max(abs(Drc_t(indOutside, :)))./max(abs(Drc_t(indMainLobe, :)))));

    powerMainLobe   = sum(abs(Drc_t(indMainLobe,:)).^2)*dr;
    powerTotal      = sum(abs(Drc_t(indTenRes,:)).^2)*dr;
    ISLR_time(kk)       = 10*log10(mean((powerTotal-powerMainLobe)./powerMainLobe));

    % The should be (almost) equal. If the bandwidth of the simulation is
    % infinite they are perfectly the same
    %figure; plot(r, db(abs(Drc_t)./max(abs(Drc_t)))); hold on; %plot(r, db(abs(Drc_f)./max(abs(Drc_f))));
    %plot([r(indMainLobe(1)), r(indMainLobe(1))], [-70, 0], 'r--');
    %plot([r(indMainLobe(end)), r(indMainLobe(end))], [-70, 0], 'r--');

    % Inversion filter (no matched filtering)
    Si      = (1./S).*rectpuls(f_ax(:)/B);

    % Make the inversion itself
    SSi     = S.*Si;

    % Back in time domain
    Drc_i   = idft(SSi,f_ax,t);

    % Normalize again the amplitudes
    Drc_i   = Drc_i*df*Nf;

    PSLR_inv(kk) = 20*log10(mean(max(abs(Drc_i(indOutside,:)))./max(abs(Drc_i(indMainLobe,:)))));

    powerMainLobe   = sum(abs(Drc_i(indMainLobe,:)).^2)*dr;
    powerTotal      = sum(abs(Drc_i(indTenRes,:)).^2)*dr;
    ISLR_inv(kk)       = 10*log10(mean((powerTotal-powerMainLobe)./powerMainLobe));
end

% Plot overlapping with the previous
% plot(r,db(abs(Drc_i)./max(abs(Drc_i))));
% xlim([-30*c/2/B 30*c/2/B]); ylim([-60 0]);
% 
% legend("Matched Filter", "Inversion Filter", 'location','northeast');
% xlabel("Range [m]");
% ylabel("Amplitude [dB]");
% title({"Range compression", sprintf("QAM with %d symbols", M)});


figure; plot(u_vect, ISLR_time, 'b-d'); hold on; 
plot(u_vect, ISLR_inv, 'b-o'); grid on; xlabel("Numerology [u]"); ylabel("ISLR [dB]");
ylim([-14 -11]);
yyaxis right
plot(u_vect, PSLR_time, 'r--d'); hold on; plot(u_vect, PSLR_inv, 'r--o'); ylabel("PSLR [dB]"); ylim([-18 -16.6]);
legend("ISLR-MF", "ISLR-ZF", "PSLR-MF", "PSLR-ZF", 'Location','best');
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
xticks(u_vect);
title("PSLR and ISLR - QPSK");
exportgraphics(gca, "C:\Users\manzoni\OneDrive - Politecnico di Milano\Conferenze\[2024] RadarConf\PSLR_ISLR_QPSK.pdf");


