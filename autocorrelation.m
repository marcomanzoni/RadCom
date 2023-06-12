function y = autocorrelation(s1,s2, dt)
%AUTOCORRELATION Summary of this function goes here
%   Detailed explanation goes here
y = conv(s1,flip(conj(s2)),"same")*dt;
end

