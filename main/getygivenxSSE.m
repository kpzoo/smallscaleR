% Risk of y infections given x introductions
function pyx = getygivenxSSE(ydom, n, pinfx)

% Assumptions and notes
% - no weighting by probability of introductions x
% - assume at least 1 infected and 1 susceptible
% - allows for super-spreading in the offspring distribution

% Matrix of y >= 0 new infections from x introductions
pyx = zeros(n+1, n+1); 

% For every number of introductions x = i-1, get every y
for i = 1:n+1
    % Resulting observed new infection probabilities
    pyx(i, :) = binopdf(ydom, n-i+1, pinfx(i));
end
