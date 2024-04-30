clear all;
close all;
w=10;
A=1;
t=-10:0.01:10
xt=A*rectpuls(t,w)
subplot(2,2,1)
plot(t,xt)
xlabel('Time')
ylabel('Amplitude')
title('Rectangular pulse')
w=-8*pi:0.01:8*pi;
for i=1:length(w)
    xw(i)=trapz(t,xt.*exp(-j*w(i).*t));
end
subplot(2,2,2)
plot(w,xw)
title('Fourier transform of rect pulse: sampling signal')
xlabel('frequency')
ylabel('Amplitude')

subplot(2,2,3)
plot(w,abs(xw))

title('magnitude spectrum')
xlabel('frequency')
ylabel('Amplitude')

subplot(2,2,4)
plot(w,angle(xw))
title('phase spectrum')
xlabel('frequency')
ylabel('Amplitude')