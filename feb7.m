clear all;
close all;
clc;
w=10;
A=1;
t=-10:0.01:10;
t0=exp(j*10*t);
xt=A*rectpuls(t,w).*t0
subplot(2,1,1)
plot(t,xt)
xlabel('Time')
ylabel('Amplitude')
title('xt')
w=-8*pi:0.01:8*pi;
for i=1:length(w)
    xw(i)=trapz(t,xt.*exp(-j*w(i).*t));
end
subplot(2,1,2)
plot(w,xw)
title('Fourier transform of xt')
xlabel('frequency')
ylabel('Amplitude')