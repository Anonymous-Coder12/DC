%% Experiment 2: Effect of sampling and verification of sampling theorem  
% Name of student:               Roll No. 33
% DOB: 27                        MOB: 11
clear 
close all
clc
 
%% Problem
Vm1=20;                     % Voltage magnitude for v1
Vm2=10;                     % Voltage magnitude for v2
DOB=27;                     % Date of birth
MOB=11;                     % Month for birth
f1=MOB*10;                  % Frequency for v1
f2=DOB*10;                  % Frequency for v2
w1=2*pi*f1;                 % Frequency for v1 (rad/sec)
w2=2*pi*f2;                 % Frequency for v2 (rad/sec)
t=0:0.00005:0.04;           % Time for plotting signal
v1=Vm1*sin(w1*t);           % v1 
v2=Vm2*sin(w2*t);           % v2
v=v1+v2;                    % Combination
 
figure(1)
subplot(2,2,1)
plot(t,v1)
hold on
plot(t,v2)
grid on
xlabel('Time (s)')
ylabel('v_1 and v_2')
title('Voltage waveforms')
legend('v_1','v_2')
 
subplot(2,2,2)
plot(t,v)
grid on
xlabel('Time (s)')
ylabel('v')
title('Combined waveform')
 
%% Defining sampling periods 
fm=max(f1,f2);          % Maximum frequency component
fs1=2*fm; Ts1=1/fs1;    % Just same
fs2=3*fm; Ts2=1/fs2;    % fs>2fm (satisfying sampling theorem)
fs3=1*fm; Ts3=1/fs3;    % fs<2fm (not satisfying sampling theorem)
 
%% fs>2fm (satisfying sampling theorem)
ts2=0:Ts2:0.04;
v1s2=Vm1*sin(w1*ts2); 
v2s2=Vm2*sin(w2*ts2);
vs2=v1s2+v2s2;
subplot(2,2,3)
stem(ts2,vs2)
xlabel('Time (s)')
ylabel('v')
title('Signal sampled with f_s>2f_m')
 
%% fs<2fm (not satisfying sampling theorem)
ts3=0:Ts3:0.04;
v1s3=Vm1*sin(w1*ts3); v2s3=Vm2*sin(w2*ts3);
vs3=v1s3+v2s3;
subplot(2,2,4)
stem(ts3,vs3)
xlabel('Time (s)')
ylabel('v')
title('Signal sampled with f_s<2f_m')
 
%% Comparison of signals (by plotting zoh type signals)
figure(2)
subplot(2,1,1)
plot(t,v,'b')
hold on
stairs(ts2,vs2,'k')
xlabel('Time (s)')
ylabel('v')
title('Comparison of signals')
legend('Original Signal','f_s>2f_m')
 
subplot(2,1,2)
plot(t,v,'b')
hold on
stairs(ts3,vs3,'--r')
xlabel('Time (s)')
ylabel('v')
title('Comparison of signals')
legend('Original Signal', 'f_s<2f_m')