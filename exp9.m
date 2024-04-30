clc;
clear all;
close all;
% Set the arrival rate (l) and service rate (u) parameters
l = 0.5;
u = 1;

% Calculate the traffic intensity (rho)
rho = l / u;

% Set the range of traffic intensity values to plot
rho_vals = linspace(0, 1, 100);

% Calculate the average delay per packet for each traffic intensity value
avg_delay_vals = zeros(size(rho_vals));
for i = 1:length(rho_vals)
    if rho_vals(i) < 1
        avg_delay_vals(i) = 1 / (u - l * rho_vals(i));
    else
        avg_delay_vals(i) = Inf;
    end
end

% Plot the results
plot(rho_vals, avg_delay_vals, 'b-', 'LineWidth', 2);
xlabel('Traffic Intensity (\rho)');
ylabel('Average Delay per Packet (seconds)');
title('M/M/1 Queue with Infinite Buffer Space');
grid on;