% Event reproduction number for a given x
function [Rx, pinfx, Vx, ysamp, Rsamp] = getRxSSEsamp(n, R0, tau, d, k, nSamps)

% Assumptions and notes
% - no weighting by probability of introductions x
% - assume at least 1 infected and 1 susceptible
% - allows for super-spreading in the offspring distribution
% - sample R0 values to obtain variance
% - get sample variations in infections y for a given x

% Reproduction numbers from x introductions
Rx = zeros(1, n+1); pinfx = Rx; Vx = Rx; 
ysamp = zeros(nSamps, n+1); Rsamp = ysamp;

% For every number of introductions x get expected infections
for i = 1:n+1
    % Sample sums of imported R0
    R0samp = gamrnd(k*(i-1), R0/k, [1 nSamps]);
    % Tranmission probability from i-1 imports for samples
    ptransRx = 1 - exp(-R0samp*tau/(d*n));
    
    % Resulting observed new infection probabilities
    if i == 1 || i == n+1
        % No imports or no susceptibles so R set to 0
         pinfx(i) = 0; Rx(i) = 0; Vx(i) = 0; 
         ysamp(:, i) = 0; Rsamp(:, i) = 0;
    else
        % Get average over the samples
        pinfx(i) = mean(ptransRx);
        % Event R as a function of x imports
        Rx(i) = ((n-i+1)/(i-1))*pinfx(i);
        % Variance around Rx from samples
        Vx(i) =  var(((n-i+1)/(i-1))*ptransRx);

        % Samples of average infections given R samples
        ysamp(:, i) = (n-i+1)*ptransRx;
        % Samples of possible R event values
        Rsamp(:, i) = ((n-i+1)/(i-1))*ptransRx;
    end


end