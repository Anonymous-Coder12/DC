% clear all
% close all
% % Creating modulating signal
% fm=2; % message signal frequency (Hz)
% n=50; % factor of sampling frequency
% K=1000;
% Ts=1/(n*fm); % sampling time
% t=0:Ts:100-Ts; % time range
% N=size(t,2);
% Fs=1/Ts; % sampling frequency
% dFs=Fs/N;
% f=-Fs/2:dFs:Fs/2-dFs;
% m=2*cos(2*pi*fm*t); % message signal
% subplot(5,1,1);
% plot(t,m);
% xlabel('Time(in s)');
% title('Modulating signal');
% % Frequency Domain
% M=fftshift(fft(m)); % FFT of the message signal
% subplot(5,1,2)
% plot(f,abs(M)/N);
% xlabel('Frequency(in hertz)');
% title('Magnitude response');
% %sound(X)
% % Pulse train generation
% T=0.2; %sampling interval
% F=1/T; % sapling frequency
% h=zeros(1,length(t)); % initialize all values to zero
% for k=-K:1:K
% h=h+(1/T)*cos(2*pi*k*F*t);
% end
% h_a=T*h/(2*K+1); % scaling 
% subplot(5,1,3);
% plot(t,h_a);
% xlabel('Time(s)');
% title('Pulse train');
% % sampling
% m_samp=m.*h_a;
% N_samp=length(m_samp);
% subplot(5,1,4); 
% plot(t,m_samp);% add the plot function 
% xlabel('Time(in s)');
% title('Sampled value');
% % Magnitude response of sampled signal
% M_samp=fftshift(fft(m_samp)); 
% subplot(5,1,5)
% plot(f,abs(M_samp)/N_samp);
% xlabel('Frequency(in hertz)');
% title('Magnitude response')

%2ndcode
% clear; clc;
% 
% f = 30e3;
% fs = 2 * f;
% Ts = 1/fs;
% t1 = 0:1e-7:5/f;
% x1 = cos(2 * pi * f * t1);
% 
% t2 = 0:Ts:5/f;
% x2 = cos(2 * pi * f * t2);
% 
% xr = zeros(size(t1));
% for i = 1:length(t1)
%     for j = 1:length(x2)
%         xr(i) = xr(i) + x2(j) * sinc(2 * fs * t1(i) - j);
%     end
% end
% 
% % Plot the original continuous signal
% subplot(3, 1, 1);
% plot(t1, x1);
% title('Continuous Signal x(t)');
% xlabel('Time');
% ylabel('Amplitude');
% 
% % Plot the sampled signal
% subplot(3, 1, 2);
% stem(t2, x2, 'r');
% title('Sampled Signal x(nT)');
% xlabel('Time');
% ylabel('Amplitude');
% 
% % Plot the reconstructed signal
% subplot(3, 1, 3);
% plot(t1, xr);
% title('Reconstructed Signal x_r(t)');
% xlabel('Time');
% ylabel('Amplitude');


%mam code
% Sampling Theorem and Signal Reconstruction
% clc;
% clear all;
% close all;
% % Parameters
% originalSamplingRate = 1000;  % Original signal sampling rate (Hz)
% minimumNyquistRate = 2 * originalSamplingRate;  % Minimum Nyquist sampling rate
% maximumSamplingRate = 10 * originalSamplingRate; % Maximum sampling rate
% 
% % Signal parameters
% originalFrequency = 50;  % Original signal frequency (Hz)
% timeVectorOriginal = 0:1/originalSamplingRate:1;  % Time vector for original signal
% 
% % Generate original signal
% originalSignal = sin(2*pi*originalFrequency*timeVectorOriginal);
% 
% % Plot the original signal
% figure;
% subplot(3, 1, 1);
% plot(timeVectorOriginal, originalSignal);
% title('Original Signal');
% xlabel('Time (s)');
% ylabel('Amplitude');
% 
% % Sample the signal at different rates
% for currentSamplingRate = [minimumNyquistRate, 5*originalSamplingRate, maximumSamplingRate]
%     samplingPeriod = 1/currentSamplingRate;  % Sampling period
%     numberOfSamples = round(1/samplingPeriod);  % Number of samples
% 
%     % Sample the signal
%     sampledSignal = originalSignal(1:numberOfSamples:end);
% 
%     % Time vector for sampled signal
%     timeVectorSampled = 0:samplingPeriod:(length(sampledSignal)-1)*samplingPeriod;
% 
%     % Reconstruct the signal using zero-order hold
%     reconstructedSignal = zeros(1, length(timeVectorOriginal));
%     reconstructedSignal(1:numberOfSamples:end) = sampledSignal;
% 
%     % Plot the sampled signal
%     subplot(3, 1, find(currentSamplingRate == [minimumNyquistRate, 5*originalSamplingRate, maximumSamplingRate])+1);
%     stem(timeVectorSampled, sampledSignal, 'r', 'Marker', 'o');
%     hold on;
%     plot(timeVectorOriginal, reconstructedSignal, 'b--');
%     hold off;
%     title(['Sampled Signal (Fs = ' num2str(currentSamplingRate) ' Hz)']);
%     xlabel('Time (s)');
%     ylabel('Amplitude');
%     legend('Sampled Signal', 'Reconstructed Signal');
% end


%divyansh
% Define the continuous-time signal
f = 100; % frequency of the signal in Hz
T = 1/f; % period of the signal in s
t = 0:0.0001:5*T; % time vector
x = sin(2*pi*f*t); % signal

% Plot the continuous-time signal
figure(1)
plot(t,x)
xlabel('Time (s)')
ylabel('x(t)')
title('Continuous-time signal')

% Sample the signal at different sampling frequencies
fs1 = 150; % sampling frequency 1 in Hz
fs2 = 200; % sampling frequency 2 in Hz
fs3 = 300; % sampling frequency 3 in Hz
Ts1 = 1/fs1; % sampling period 1 in s
Ts2 = 1/fs2; % sampling period 2 in s
Ts3 = 1/fs3; % sampling period 3 in s
n1 = 0:Ts1:5*T; % sample index vector 1
n2 = 0:Ts2:5*T; % sample index vector 2
n3 = 0:Ts3:5*T; % sample index vector 3
x1 = sin(2*pi*f*n1); % sampled signal 1
x2 = sin(2*pi*f*n2); % sampled signal 2
x3 = sin(2*pi*f*n3); % sampled signal 3

% Plot the sampled signals
figure(2)
subplot(3,1,1)
stem(n1,x1)
hold on
plot(n1,x1)
hold off
xlabel('Time (s)')
ylabel('x[n]')
title(['Sampled signal 1 with fs = ',num2str(fs1),' Hz'])
subplot(3,1,2)
stem(n2,x2)
hold on
plot(n2,x2)
hold off
xlabel('Time (s)')
ylabel('x[n]')
title(['Sampled signal 2 with fs = ',num2str(fs2),' Hz'])
subplot(3,1,3)
stem(n3,x3)
hold on
plot(n3,x3)
hold off
xlabel('Time (s)')
ylabel('x[n]')
title(['Sampled signal 3 with fs = ',num2str(fs3),' Hz'])

% Reconstruct the signal using sinc interpolation
xr1 = zeros(size(t)); % reconstructed signal 1
xr2 = zeros(size(t)); % reconstructed signal 2
xr3 = zeros(size(t)); % reconstructed signal 3
for i = 1:length(n1)
    xr1 = xr1 + x1(i)*sinc((t-n1(i))/Ts1); % sinc interpolation formula
end
for i = 1:length(n2)
    xr2 = xr2 + x2(i)*sinc((t-n2(i))/Ts2); % sinc interpolation formula
end
for i = 1:length(n3)
    xr3 = xr3 + x3(i)*sinc((t-n3(i))/Ts3); % sinc interpolation formula
end

% Plot the reconstructed signals
figure(3)
subplot(3,1,1)
plot(t,xr1)
xlabel('Time (s)')
ylabel('xr1(t)')
title(['Reconstructed signal 1 with fs = ',num2str(fs1),' Hz'])
subplot(3,1,2)
plot(t,xr2)
xlabel('Time (s)')
ylabel('xr2(t)')
title(['Reconstructed signal 2 with fs = ',num2str(fs2),' Hz'])
subplot(3,1,3)
plot(t,xr3)
xlabel('Time (s)')
ylabel('xr3(t)')
title(['Reconstructed signal 3 with fs = ',num2str(fs3),' Hz'])