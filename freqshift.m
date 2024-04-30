% Frequency Shifting Property
% Original signal
t = -1:0.01:1;
x = rectpuls(t, 0.5); % Rectangular pulse
% Frequency shift
f0 = 2;
x_shifted_frequency = x .* exp(1i * 2 * pi * f0 * t); % Shifted frequency
% Fourier Transform of the original signal
X = fft(x);
% Fourier Transform of the shifted frequency signal
X_shifted_frequency = fft(x_shifted_frequency);
% Plot the original and frequency-shifted signals in time domain
subplot(2, 1, 1);
plot(t, x, 'b', t, real(x_shifted_frequency), 'r--');
title('Time Domain - Original and Frequency-Shifted Signals');
legend('Original', 'Frequency-Shifted');
xlabel('Time');
ylabel('Amplitude');
% Plot the Fourier Transforms in frequency domain
subplot(2, 1, 2);
f = linspace(-50, 50, length(X));
plot(f, abs(fftshift(X)), 'b', f, abs(fftshift(X_shifted_frequency)), 'r--');
title('Frequency Domain - Fourier Transforms');
legend('Original', 'Frequency-Shifted');
xlabel('Frequency');
ylabel('Magnitude');
