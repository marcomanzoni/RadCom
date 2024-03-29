function alpha = snowPowerAttenuation(f)
%SNOWPOWERATTENUATION Taken from: 
% Complex dielectric permittivity measurements fromground-penetrating radar 
% data to estimate snowliquid water content in the pendular regime

mu0         = 1.25663706e-6;
e0          = 8.85418782e-12;
rho_d_kg_m3 = 275; % kg/m^3
rho_d       = rho_d_kg_m3/1000; % g/cm^3
W           = 0.01;
ew_prime    = 80;
ew_second   = 13;

ed = (1 + 1.7*rho_d + 0.7*rho_d^2);
es_prime = (0.1*W + 0.8*W^2)*ew_prime + ed;
es_second = (0.1*W + 0.8*W^2)*ew_second; 

alpha = sqrt(mu0/(es_prime*e0))*(es_second*e0/2)*2*pi*f;

attenuation_db_m = -alpha*20*log10(exp(1)) % db/m of loss

end

