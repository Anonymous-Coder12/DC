% Given generator matrix G
G = [1 1 0 0 1 0 1; 
     0 1 1 1 1 0 0; 
     1 1 1 0 0 1 1];

% Initialize G_prime with zeros
G_prime = zeros(size(G));

% Perform row operations to get systematic form
for i = 1:size(G, 1)
    % Find the index of the leading 1 in the current row
    leading_1_index = find(G(i,:) == 1, 1);
    
    % Update G_prime with the current row
    G_prime(i, :) = G(i, :);
    
    % Perform row operations to get zeros below the leading 1
    for j = 1:size(G, 1)
        if j ~= i
            factor = G(j, leading_1_index);
            G(j, :) = mod(G(j, :) - factor * G(i, :), 2);
        end
    end
end

% Display the systematic form of the generator matrix G'
disp('Generator Matrix G'' for Systematic Code C''');
disp(G_prime);