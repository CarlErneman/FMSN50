clc;
close all;
clear;
load coal_mine_disasters.mat;


%% 1c)
% Define our interval for from the data given
t_start = 1658;
t_end = 1980;

% Number of samples (N), breakpoints (d) and burn-in (burn_in)
N = 25000;
d = 6; %Was 6
burn_in = 10000;

% Assign value to the hyperparameter (Psi) 
psi = 20;

lambda_save = zeros(d,d-1);
breakpoint_save = zeros(d-1, d);

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
        [~,t] = MH_algorithm(t, lambda, tau, rho);
    end
    
    %Run after burn-in for samlpes
    for j = 1:N
        % Draw theta from corresponding marginal posterior distribution
        theta =  gamrnd(2*length(lambda) + 2, 1./(psi + sum(lambda)));
        % Draw lambda from corresponding marginal posterior distribution
        lambda = posterior_lambda(theta, t, tau);
        lambda_save(1:i, i-1) = lambda;
        [accepted_proposals,t] = MH_algorithm(t, lambda, tau, rho);
        breakpoints_update(j, :) = t;
    end
    breakpoint_save (i-1, 1:i) = breakpoints_update(25000, 1:i);

    figure;
    plot(breakpoints_update);
    xlabel("Itterations");
    ylabel("Year");
    title("Breakpoint movement for d=" + num2str(i-1));
    set(gca, 'Fontsize', 10);
    %filename = "breakpont_itter_d" + num2str(i-1) + ".png";
    %saveas(gcf,filename)

    figure;
    for k = 2:i
        histogram(breakpoints_update(:,k), 20)
        hold on
    end
    title("The histogram of the chain for d = " + num2str(k-1)+  " of breakpoints")
    xlabel('Breakpoints')
    ylabel('Intensity')
    set(gca, 'Fontsize', 10);
    %filename = "histogram_itter_d" + num2str(i-1) + ".png";
    %saveas(gcf,filename)
end


% Plot the data
figure;
hold on
plot(tau, 1:751, "b", "LineWidth", 2);
xlabel("Year");
ylabel("Number of disasters");
title("Total number of disasters during 1658-1980 and marked brakpoints using d=" + num2str(d-1));
set(gca, 'Fontsize', 10);
xline(breakpoints_update(25000, 2:d), "--s")
hold off

%filename = "disasters.png";
%saveas(gcf,filename)

%% 1d)
% Initialize as in 1c) 
load coal_mine_disasters.mat;
% Define our interval for from the data given
t_start = 1658;
t_end = 1980;

% Number of samples (N), breakpoints (d), rho and burn-in (burn_in)
N = 25000;
d = 6;
burn_in = 10000;
rho = 0.01 * ones(d,1);

% Construct the breakpoints-vector 
step_size = (t_end - t_start)/d;
t_middle = t_start:step_size:t_end;
t = [t_start,t_middle(2:end-1),t_end];
breakpoints_update = zeros(N, length(t));

psi = 50; 

%Preallocating posteriors for speed 
theta_mean = zeros(psi, 1);
theta_var = zeros(psi, 1);

lambda_mean = zeros(psi, d);
lambda_var = zeros(psi, d);

t_mean = zeros(psi, length(t));
t_var = zeros(psi, length(t));

% Implement same method as earlier, but looping over different psi's now
% instead of breakpoints
for i = 1:psi
    % Assign the distributions to hyperprior and prior, 
    %taking the inverse of the second parameter due to MATLABs notation 
    theta = gamrnd(2, 1/i); % Gamma(2, Psi)
    lambda = gamrnd(2, 1/theta, 1, d); % Gamma(2, theta)

    % Save temporary values of theta and lambda to calculate the mean &
    % variance for each value on psi
    theta_temp = zeros(N, 1);
    lambda_temp = zeros(N,d);

    % start burn-in to reach stationary distribution "sufficiently close"
    for j = 1:burn_in
        % Draw theta from corresponding marginal posterior distribution
        theta =  gamrnd(2*length(lambda) + 2, 1./(i + sum(lambda)));
        % Draw lambda from corresponding marginal posterior distribution
        lambda = posterior_lambda(theta, t, tau);
        % Discard the accepted proposals
        [~,t] = MH_algorithm(t, lambda, tau, rho);
    end
    % Same thing again but now for our samples N
    % but now assigning accepted proposals from MH-algorithm
    for j = 1:N
        % Draw theta from corresponding marginal posterior distribution
        theta =  gamrnd(2*length(lambda) + 2, 1./(i + sum(lambda)));
        % Draw lambda from corresponding marginal posterior distribution
        lambda = posterior_lambda(theta, t, tau);
        [accepted_proposals,t] = MH_algorithm(t, lambda, tau, rho);

        theta_temp(j) = theta;
        lambda_temp(j, :)= lambda;
        % Save the time for the event
        breakpoints_update(j,:) = t;
    end
    % Compute the means for psi, before incrementing
    theta_mean(i) = mean(theta_temp);
    theta_var(i) = var(theta_temp);

    lambda_mean(i, :) = mean(lambda_temp);
    lambda_var(i, :) = var(lambda_temp);

    t_mean(i, :) = mean(breakpoints_update);
    t_var(i, :) = var(breakpoints_update);
end
%% Plots 1d)
figure();
plot(theta_mean, "LineWidth", 2, 'Marker','+');
xlabel("Different values on \Psi");
ylabel("Mean"); 
title("\theta's dependence on different values of \Psi by comparing mean");
grid on

figure();
plot(theta_var,"LineWidth", 2, 'Marker','x');
xlabel("Different values on \Psi");
ylabel("Variance"); 
title("\theta's dependence on different values of \Psi by comparing variance");
grid on

figure();
plot(lambda_mean,"LineWidth", 2, 'Marker','x');
xlabel("Different values on \Psi");
ylabel("Mean"); 
title("\lambda's dependence on different values of \Psi by comparing mean");
legend(['\lambda_1'; '\lambda_2'; '\lambda_3'; '\lambda_4'; '\lambda_5'])
%grid on

figure();
plot(lambda_var,"LineWidth", 2, 'Marker','x');
xlabel("Different values on \Psi");
ylabel("Variance"); 
title("\lambda's dependence on different values of \Psi by comparing variance");
legend(['\lambda_1'; '\lambda_2'; '\lambda_3'; '\lambda_4'; '\lambda_5'])
%grid on

figure();
plot(t_mean(:,2:end-1), 'Marker','x');
xlabel("Different values on \Psi");
ylabel("Mean"); 
title("Breakpoints dependence on different values of \Psi by comparing mean");
legend(['t_1'; 't_2'; 't_3'; 't_4'; 't_5'])

figure();
plot(t_var(:,2:end-1), 'Marker','x');
xlabel("Different values on \Psi");
ylabel("Variance"); 
title("Breakpoints dependence on different values of \Psi by comparing variance");
legend(['t_1'; 't_2'; 't_3'; 't_4'; 't_5'])

%% 1e)
psi = 20;
% Range for the values of rho tested
rhos = 0.01:0.001:0.1;
% Pre-allocate for speed
rho = zeros(d, length(rhos));
accepted_proposals_rhos = zeros(length(rhos), d-1);
% Assign the values from rhos to the matrix
for i = 1:d 
    rho(i,:) = rhos; 
end

for i = 1:length(rhos)
    % Assign the distributions to hyperprior and prior, 
    %taking the inverse of the second parameter due to MATLABs notation 
    theta = gamrnd(2, 1/psi); % Gamma(2, Psi)
    lambda = gamrnd(2, 1/theta, 1, d); % Gamma(2, theta)
    
    % Construct the breakpoints-vector 
    step_size = (t_end - t_start)/d;
    t_middle = t_start:step_size:t_end;
    t = [t_start,t_middle(2:end-1),t_end];
    breakpoints_update = zeros(N, length(t));
    
    % start burn-in to reach stationary distribution "sufficiently close"
    for j = 1:burn_in
        % Draw theta from corresponding marginal posterior distribution
        theta =  gamrnd(2*length(lambda) + 2, 1./(psi + sum(lambda)));
        % Draw lambda from corresponding marginal posterior distribution
        lambda = posterior_lambda(theta, t, tau);
        % Discard the accepted proposals
        [~,t] = MH_algorithm(t, lambda, tau, rho(1,i));
    end
    
    % Same thing again but now for our samples N
    % but now assigning accepted proposals from MH-algorithm
    for j = 1:N
        % Draw theta from corresponding marginal posterior distribution
        theta =  gamrnd(2*length(lambda) + 2, 1./(psi + sum(lambda)));
        % Draw lambda from corresponding marginal posterior distribution
        lambda = posterior_lambda(theta, t, tau);
        [accepted_proposals,t] = MH_algorithm(t, lambda, tau, rho(1,i));
        % Save accepted proposals for different rhos
        accepted_proposals_rhos(i,:) = accepted_proposals_rhos(i,:) + accepted_proposals;
    end
end
%%
% Sum all the rows to find the rate
accepted_proposals_rate = sum(accepted_proposals_rhos,2) / (N*(d-1));

figure();
plot(rhos, accepted_proposals_rate, "LineWidth",1.5, "Color", "black");
yline(0.3, "LineWidth", 1, "Color","yellow");
yline([0.23 0.44],"--g","LineWidth", 1);
ylabel("Acceptance rate");
xlabel("\rho");
xlim([0.01 0.1]);
title("Acceptance rate for different values of \rho");
legend("Acceptance rate curve", "Optimal acceptance rate",...
    "Rule of thimb range for acceptance rate");

%% Simulation for autocorrelation 
% Check for four different values of rhos, found on previous plot
rhos = [.01, .02, .03, .04];
rho = zeros(d, length(rhos));
% Assign the values from rhos to the matrix
for i = 1:d 
    rho(i,:) = rhos; 
end

for i = 1:length(rhos)
    % Assign the distributions to hyperprior and prior, 
    %taking the inverse of the second parameter due to MATLABs notation 
    theta = gamrnd(2, 1/psi); % Gamma(2, Psi)
    lambda = gamrnd(2, 1/theta, 1, d); % Gamma(2, theta)
    
    % Construct the breakpoints-vector 
    step_size = (t_end - t_start)/d;
    t_middle = t_start:step_size:t_end;
    t = [t_start,t_middle(2:end-1),t_end];
    breakpoints_update = zeros(N, length(t));
    
    % start burn-in to reach stationary distribution "sufficiently close"
    for j = 1:burn_in
        % Draw theta from corresponding marginal posterior distribution
        theta =  gamrnd(2*length(lambda) + 2, 1./(psi + sum(lambda)));
        % Draw lambda from corresponding marginal posterior distribution
        lambda = posterior_lambda(theta, t, tau);
        % Discard the accepted proposals
        [~,t] = MH_algorithm(t, lambda, tau, rho(1,i));
    end
    
    % Same thing again but now for our samples N
    % but now assigning accepted proposals from MH-algorithm
    for j = 1:N
        % Draw theta from corresponding marginal posterior distribution
        theta =  gamrnd(2*length(lambda) + 2, 1./(psi + sum(lambda)));
        % Draw lambda from corresponding marginal posterior distribution
        lambda = posterior_lambda(theta, t, tau);
        [accepted_proposals,t] = MH_algorithm(t, lambda, tau, rho(1,i));
        breakpoints_update(j,:) = t;
    end
    % Plot the autocorrelation function for the different rhos
    figure();
    subplot(2,3,1); % t_1
    autocorr(breakpoints_update(:,2),5e2);

    subplot(2,3,2); % t_2
    autocorr(breakpoints_update(:,3),5e2);

    subplot(2,3,3); % t_3
    autocorr(breakpoints_update(:,4),5e2);

    subplot(2,3,4); % t_4
    autocorr(breakpoints_update(:,5),5e2);

    subplot(2,3,5); % t_5
    autocorr(breakpoints_update(:,6),5e2);
end

%% 2b)

load("atlantic.txt");

f_inv = @(u, mu, beta) mu - beta*log(-log(u));
n = length(atlantic);

% Choose how many new data sets we want to generate
B = 1000;

% Estimate the parameters mu and beta with the est_gumble function
[beta_hat, mu_hat] = est_gumbel(atlantic);

beta_boot = zeros(1,B);
mu_boot = zeros(1,B);

% Perform the bootstrap method
for i = 1:B
    y_boot = f_inv(rand(n,1), mu_hat, beta_hat);
    [beta_temp, mu_temp] = est_gumbel(y_boot);
    beta_boot(i) = beta_temp;
    mu_boot(i) = mu_temp;
end

% sorting to obtain quantiles for beta and mu
beta_delta = sort(beta_boot-beta_hat);
mu_delta = sort(mu_boot-mu_hat);

% CB level
alpha = 0.01;

% Construct the two confidence intervals 
beta_lower = beta_hat-beta_delta(ceil((1-alpha/2)*B)); 
beta_upper = beta_hat+beta_delta(ceil((1-alpha/2)*B));

mu_lower = mu_hat-mu_delta(ceil((1-alpha/2)*B)); 
mu_upper = mu_hat+mu_delta(ceil((1-alpha/2)*B));

