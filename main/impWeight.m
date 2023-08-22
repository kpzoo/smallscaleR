% Weight the samples of event R or y infections by P(x)
function [yw, Rw] = impWeight(n, ysamp, Rsamp, Px, nSamps)

% All y and R samples reweighted
yw = cell(1, n+1); Rw = yw;

% Reweight the samples with the import probability
for i = 1:n+1
    % Bernoulli sample according to Px
    ids = binornd(1, Px(i), [1 nSamps]);
    ids = find(ids);
    % Select column for this x and thin
    yw{i} = ysamp(:, i); Rw{i} = Rsamp(:, i);
    yw{i} = yw{i}(ids); Rw{i} = Rw{i}(ids);
end

% Form a distribution to marginalise x
yw = cell2mat(yw')'; Rw = cell2mat(Rw')';