function alpha = snowPowerAttenuationLuca(f,W)
% SNOWPOWERATTENUATION Taken from: 
% Complex dielectric permittivity measurements fromground-penetrating radar 
% data to estimate snowliquid water content in the pendular regime
%
% Model:
%   es_prime    = (0.1*W + 0.8*W^2)*ew_prime + ed
%   es_second   = (0.1*W + 0.8*W^2)*ew_second

mu0         = 1.25663706e-6;
e0          = 8.85418782e-12;
rho_d_kg_m3 = 500;                  % kg/m^3
rho_d       = rho_d_kg_m3/1000;     % g/cm^3
ed          = (1 + 1.7*rho_d + 0.7*rho_d^2);


% Water parameters
% "The Complex Dielectric Constant of Snow at Microwave Frequencies"
ew_inf      = 4.9;      % Infinite frequency dielectric constant at 0 °C
ew_s        = 87.74;    % Static dielectric constant at 0 °C
tau_w       = 18e-12;   % Relaxation time -> [17,18] ps

% Debye relaxation equation -> ew as a function of frequency
ew          = ew_inf + (ew_s-ew_inf)/(1+1j*2*pi*f*tau_w);
ew_prime    = real(ew);
ew_second   = -imag(ew);

if 1
    % Method from
    % "A Microwave Frequency Range Experiment for the Measurement of Snow 
    % Density and Liquid Water Content"
    es_prime    = (0.1*W + 0.8*W^2)*ew_prime + (1 + 1.7*(rho_d-W) + 0.7*(rho_d-W)^2);
    es_second   = (0.1*W + 0.8*W^2)*ew_second;

else
    % Mixing parameters from 
    % "Microwave Dielectric Properties of Surface Snow"
    chi         = 1/4;      % Coupling factor           -> need more data!!
    A0          = .003;     % Depolarization factor     -> need more data!!

    % Mixing equation
    es          = ed + chi*W*(ew-ed)/(1+(ew/ed-1)*A0)
    es_prime    = real(es);
    es_second   = -imag(es);
end


alpha = sqrt(mu0/(es_prime*e0))*(es_second*e0/2)*2*pi*f;

attenuation_db_m = -alpha*20*log10(exp(1)); % db/m of loss
fprintf("Loss: %.3f dB/m\n", attenuation_db_m);

end

