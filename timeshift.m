% Time Shifting Property
% Original signal
t = -1:0.01:1;
x = rectpuls(t, 0.5); % Rectangular pulse
% Time shift
t0 = 0.5;
x_shifted = rectpuls(t - t0, 0.5); % Shifted pulse
% Fourier Transform of the original signal
X = fft(x);
% Fourier Transform of the shifted signal
X_shifted = fft(x_shifted);
% Plot the original and shifted signals in time domain
subplot(2, 1, 1);
plot(t, x, 'b', t, x_shifted, 'r--');
title('Time Domain - Original and Shifted Signals');
legend('Original', 'Shifted');
xlabel('Time');
ylabel('Amplitude');
% Plot the Fourier Transforms in frequency domain
subplot(2, 1, 2);
f = linspace(-50, 50, length(X));
plot(f, abs(fftshift(X)), 'b', f, abs(fftshift(X_shifted)), 'r--');
title('Frequency Domain - Fourier Transforms');
legend('Original', 'Shifted');
xlabel('Frequency');
ylabel('Magnitude');