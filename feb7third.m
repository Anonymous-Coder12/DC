clc;
clear all;
close all;
uniform_data = rand(1, 1000);

mu = 0;
sigma = 1;
gaussian_data = mu + sigma * randn(1, 1000);

figure;
subplot(2, 2, 1);
histogram(uniform_data, 'Normalization', 'pdf');
title('Uniform Random Number - PDF');
xlabel('Value');
ylabel('Probability');

subplot(2, 2, 2);
histogram(uniform_data, 'Normalization', 'cdf');
title('Uniform Random Number - CDF');
xlabel('Value');
ylabel('Cumulative Probability');

subplot(2, 2, 3);
histogram(gaussian_data, 'Normalization', 'pdf');
title('Gaussian Random Number - PDF');
xlabel('Value');
ylabel('Probability');

subplot(2, 2, 4);
histogram(gaussian_data, 'Normalization', 'cdf');
title('Gaussian Random Number - CDF');
xlabel('Value');
ylabel('Cumulative Probability');

mean_uniform = mean(uniform_data);
variance_uniform = var(uniform_data);

mean_gaussian = mean(gaussian_data);
variance_gaussian = var(gaussian_data);

disp('Uniform Random Numbers:');
disp(['Mean: ' num2str(mean_uniform)]);
disp(['Variance: ' num2str(variance_uniform)]);

disp('Gaussian Random Numbers:');
disp(['Mean: ' num2str(mean_gaussian)]);
disp(['Variance: ' num2str(variance_gaussian)]);