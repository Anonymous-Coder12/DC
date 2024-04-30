clc;
clear all;
close all;
% Stop-and-Wait Protocol Simulation

% Parameters
num_packets = 10; % Number of packets to be transmitted
p_error = 0.1;    % Probability of packet error

% Simulation
successful_transmissions = 0;
total_transmissions = 0;

for i = 1:num_packets
    % Sender transmits packet
    fprintf('Sender: Packet %d transmitted.\n', i);
    
    % Simulate channel
    if rand > p_error
        % Packet successfully received
        fprintf('Receiver: Packet %d received successfully.\n', i);
        successful_transmissions = successful_transmissions + 1;
    else
        % Packet lost or corrupted
        fprintf('Receiver: Packet %d lost or corrupted.\n', i);
    end
    
    % Sender waits for acknowledgement
    fprintf('Sender: Waiting for acknowledgement...\n');
    pause(1); % Simulating delay
    
    % Simulate acknowledgement
    if rand > p_error
        % Acknowledgement received
        fprintf('Receiver: Acknowledgement for Packet %d received.\n', i);
        total_transmissions = total_transmissions + 1;
    else
        % Acknowledgement lost or corrupted
        fprintf('Receiver: Acknowledgement for Packet %d lost or corrupted.\n', i);
    end
    
    % Sender proceeds to next packet
    fprintf('\n');
end

% Performance Evaluation
fprintf('Performance Evaluation:\n');
fprintf('Successful Transmissions: %d\n', successful_transmissions);
fprintf('Total Transmissions: %d\n', total_transmissions);
fprintf('Packet Error Rate: %.2f\n', (num_packets - successful_transmissions) / num_packets);