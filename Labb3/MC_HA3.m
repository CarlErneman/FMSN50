clc;
close all;
clear;
load coal_mine_disasters.mat;

% Plot the data
plot(tau, 1:751, "b", "LineWidth", 2);
xlabel("Year");
ylabel("Number of disasters");
title("Total number of disasters during 1658-1980");
figure(1);

%% 1c)
% Define our interval for from the data given
t_start = 1658;
t_end = 1980;

% Number of samples (N), breakpoints (d) and burn-in (burn_in)
N = 25000;
d = 6;
burn_in = 10000;

% Assign value to the hyperparameter (Psi) 
psi = 20;

for i = 2:d
    % Assign the distributions to hyperprior and prior, 
    %taking the inverse of the second parameter due to MATLABs notation 
    theta = gamrnd(2, 1/psi); % Gamma(2, Psi)
    lambda = gamrnd(2, 1/theta, 1, i); % Gamma(2, theta)
    rho = 0.01 * ones(1,d);
    
    % Construct the breakpoints-vector 
    step_size = (t_end - t_start)/i;
    t_middle = t_start:step_size:t_end;
    t = [t_start,t_middle(2:end-1),t_end];
    breakpoints_update = zeros(N, length(t));
    
    % start burn-in to reach stationary distribution "sufficiently close"
    for j = 1:burn_in
        % Draw theta from corresponding marginal posterior distribution
        theta =  gamrnd(2*length(lambda) + 2, 1./(psi + sum(lambda)));
        % Draw lambda from corresponding marginal posterior distribution
        lambda = posterior_lambda(theta, t, tau);
        [accepted_proposals,t] = MH_algorithm(t, lambda, tau, rho);
    end
end