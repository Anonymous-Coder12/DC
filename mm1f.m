clc;
clear all;
close all;
% Set system parameters
lambda = 0.5; % arrival rate
mu = 1; % service rate

% Define buffer space values
u_values = 1:10;

% Calculate and plot blocking probability for each buffer space value
figure;
hold on;
for u = u_values
    % Compute traffic intensity and blocking probability
    rho = lambda / mu;
    p0 = 1 / (1 + sum((rho .^ (0:u)) ./ factorial(0:u)));
    p_blocking = (rho .^ u / factorial(u)) * p0;
    
    % Plot the blocking probability for this buffer space value
    plot(u, p_blocking, 'bo');
end

% Set axis labels and title
xlabel('Buffer Space (u)');
ylabel('Blocking Probability');
title('M/M/1 Queue with Finite Buffer Space');

% Add a line to show the theoretical maximum for comparison
rho = lambda / mu;
p_blocking_max = rho ^ (u_values(end) + 1);
plot(u_values, p_blocking_max * ones(size(u_values)), 'k--');

% Add a legend
legend('Blocking Probability', 'Theoretical Maximum', 'Location', 'NorthWest');

hold off;
