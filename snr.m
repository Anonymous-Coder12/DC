clc;
clear all;
close all;
% Set the parameters
N = 1000; % number of samples
fs = 1000; % sampling frequency (Hz)
f = 10; % signal frequency (Hz)
A = 1; % signal amplitude
Q = 10; % maximum quantization level

% Generate the signal
t = (0:N-1)'/fs; % time vector
x = A*sin(2*pi*f*t); % sinusoidal signal

% Compute the SQNR for each quantization level
sqnr = zeros(Q, 1);
for q = 1:Q
    % Quantize the signal
    xq = round(x*(2^(q-1))/A)*(A/(2^(q-1)));
    disp(xq);
    % Compute the quantization error and the SQNR
    e = x - xq;
    sqnr(q) = 10*log10(var(x)/var(e));
end

% Plot the SQNR versus quantization levels
plot(1:Q, sqnr, 'o-')
xlabel('Quantization Levels')
ylabel('Signal-to-Quantization-Noise Ratio (dB)')
title('Uniform Quantization SQNR vs. Quantization Levels')