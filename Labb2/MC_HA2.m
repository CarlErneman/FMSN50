clc; clear;



%% 3). Naive approach

disp("3). Sequential Importance Sampling, Naive approach")

n = 10 + 1; %Nbr of steps +1(to ignore the initial state)
d = 2; %Dimensions
N = 10000; %Nbr of particles


X = zeros(n, d, N);
dir_mat = [eye(d); -1*eye(d)];
N_sa = zeros(n,1);
c_n3 = zeros(n,1);


for particle = 1:N
    for stepnbr = 2:n
        nextX = datasample(dir_mat,1);
        X(stepnbr, :, particle) = X(stepnbr-1, :, particle) + nextX;

        %Check new step against all old steps. If not SAW, increase
        %error counter
        if ~ismember(X(stepnbr, :, particle), X(1:stepnbr-1, :, particle), 'rows')
            N_sa(stepnbr, 1) = N_sa(stepnbr, 1) + 1;
        else
            break;
        end
    end
end

for stepnbr = 2:n
    c_n3(stepnbr,1) = (2*d)^(stepnbr-1)*N_sa(stepnbr, 1)/N;
end

fprintf('N_sa for %i steps: \r\n', n-1)
disp(N_sa(2:n, 1))

fprintf('c_n from 1 to %i steps: \r\n', n-1)
disp(c_n3(2:n, 1))




%% 4). 

disp("___________________________________________________________")
disp("Improved Sequential Importance Sampling, g=SAW")

n = 50 + 1; %Nbr of steps +1(to ignore the initial state)
d = 2; %Dimensions
N = 10000; %Nbr of particles


X = zeros(n, d, N);
dir_mat = [eye(d); -1*eye(d)];
N_sa = zeros(n,1);
c_n4 = zeros(n,1);
w_i = [ones(1, N); zeros(n, N)];


for particle = 1:N
    for stepnbr = 1:n
        X_ki_steps = X(stepnbr, :, particle) + dir_mat;
        X_0ki = X(1:stepnbr, :, particle);
        free_coordinates = setdiff(X_ki_steps, X_0ki, 'rows');

        if isempty(free_coordinates)
            X(stepnbr+1, :, particle) = X(stepnbr, :, particle);
            w_i(stepnbr+1, particle) = 0;
        else
            nextX = datasample(free_coordinates,1);
            X(stepnbr+1, :, particle) = nextX;
            nextXnbr = size(free_coordinates, 1);
            w_i(stepnbr+1, particle) = w_i(stepnbr, particle)*nextXnbr;
        end
    end
end

for stepnbr = 2:n
    c_n4(stepnbr,1) = mean(w_i(stepnbr, :));
end

fprintf('c_n from 1 to %i steps: \r\n', n-1)
disp(c_n4(2:n, 1))

figure(1)
h1 = histogram(w_i(:,1)./sum(w_i(:,1)));
set(gca, 'Xscale',  'log')

%% 5). Sequential Importance Sampling With Resampling

disp("___________________________________________________________")
disp("Sequential Importance Sampling With Resampling")


n = 10 + 1; %Nbr of steps +1(to ignore the initial state)
d = 2; %Dimensions
N = 10000; %Nbr of particles


X = zeros(n, d, N);
dir_mat = [eye(d); -1*eye(d)];
N_sa = zeros(n,1);

c_n5 = zeros(n,1);
c_n5(1) = 1;
w_i = ones(n, N);


for stepnbr = 1:n
    ind = randsample(N, N, true, w_i(stepnbr, :));
    X(1:stepnbr, :, :) = X(1:stepnbr, :, ind);
    for particle = 1:N
        X_ki_steps = X(stepnbr, :, particle) + dir_mat;
        X_0ki = X(1:stepnbr, :, particle);
        free_coordinates = setdiff(X_ki_steps, X_0ki, 'rows');

        if isempty(free_coordinates)
            X(stepnbr+1, :, particle) = X(stepnbr, :, particle);
            w_i(stepnbr+1, particle) = 0;
        else
            nextX = datasample(free_coordinates,1);
            X(stepnbr+1, :, particle) = nextX;
            nextXnbr = size(free_coordinates, 1);
            w_i(stepnbr+1, particle) = nextXnbr;
        end
    end
end

for stepnbr = 1:n
    c_n5(stepnbr+1,1) = c_n5(stepnbr,1)*mean(w_i(stepnbr+1, :));
end

%c_n5 = round(c_n5)
%c_n5

%% 6).

disp("___________________________________________________________")
disp("Estimating Parameters")

replicates = 10;
A = zeros(replicates, 1);
mu = zeros(replicates, 1);
gamma = zeros(replicates, 1);

n = 100 + 1; %Nbr of steps +1(to ignore the initial state)
d = 2; %Dimensions
N = 10000; %Nbr of particles

for i = 1:replicates

    X = zeros(n, d, N);
    dir_mat = [eye(d); -1*eye(d)];
    N_sa = zeros(n,1);
    
    c_n6 = zeros(n,1);
    c_n6(1) = 1;
    w_i = ones(n, N);
    
    
    for stepnbr = 1:n
        ind = randsample(N, N, true, w_i(stepnbr, :));
        X(1:stepnbr, :, :) = X(1:stepnbr, :, ind);
        for particle = 1:N
            X_ki_steps = X(stepnbr, :, particle) + dir_mat;
            X_0ki = X(1:stepnbr, :, particle);
            free_coordinates = setdiff(X_ki_steps, X_0ki, 'rows');
    
            if isempty(free_coordinates)
                X(stepnbr+1, :, particle) = X(stepnbr, :, particle);
                w_i(stepnbr+1, particle) = 0;
            else
                nextX = datasample(free_coordinates,1);
                X(stepnbr+1, :, particle) = nextX;
                nextXnbr = size(free_coordinates, 1);
                w_i(stepnbr+1, particle) = nextXnbr;
            end
        end
    end
    for stepnbr = 1:n
        c_n6(stepnbr+1,1) = c_n6(stepnbr,1)*mean(w_i(stepnbr+1, :));
    end

    step_vect = 1:1:n;
    cn_log = log(c_n6(step_vect+1));
    x_1 = ones(length(step_vect), 1);
    x_2 = step_vect';
    x_3 = log(step_vect)';
    beta = regress(cn_log, [x_1, x_2, x_3]);
    A(i) = exp(beta(1));
    mu(i) = exp(beta(2));
    gamma(i) = beta(3) + 1;
end

A_var = var(A);
mu_var = var(mu);
gamma_var = var(gamma);

%% 9).

replicates = 10;
A5 = zeros(replicates, 1);
mu5 = zeros(replicates, 1);
gamma5 = zeros(replicates, 1);

n = 100 + 1; %Nbr of steps +1(to ignore the initial state)
d = 5; %Dimensions
N = 10000; %Nbr of particles

for i = 1:replicates
    X = zeros(n, d, N);
    dir_mat = [eye(d); -1*eye(d)];
    N_sa = zeros(n,1);
    
    c_n9 = zeros(n,1);
    c_n9(1) = 1;
    w_i = ones(n, N);
    
    
    for stepnbr = 1:n
        ind = randsample(N, N, true, w_i(stepnbr, :));
        X(1:stepnbr, :, :) = X(1:stepnbr, :, ind);
        for particle = 1:N
            X_ki_steps = X(stepnbr, :, particle) + dir_mat;
            X_0ki = X(1:stepnbr, :, particle);
            free_coordinates = setdiff(X_ki_steps, X_0ki, 'rows');
    
            if isempty(free_coordinates)
                X(stepnbr+1, :, particle) = X(stepnbr, :, particle);
                w_i(stepnbr+1, particle) = 0;
            else
                nextX = datasample(free_coordinates,1);
                X(stepnbr+1, :, particle) = nextX;
                nextXnbr = size(free_coordinates, 1);
                w_i(stepnbr+1, particle) = nextXnbr;
            end
        end
    end
    for stepnbr = 1:n
        c_n9(stepnbr+1,1) = c_n9(stepnbr,1)*mean(w_i(stepnbr+1, :));
    end

    step_vect = 1:1:n;
    cn_log = log(c_n9(step_vect+1));
    x_1 = ones(length(step_vect), 1);
    x_2 = step_vect';
    x_3 = log(step_vect)';
    beta = regress(cn_log, [x_1, x_2, x_3]);
    A5(i) = exp(beta(1));
    mu5(i) = exp(beta(2));
    gamma5(i) = beta(3) + 1;
end

mu5_asymp = 2*d-1-(1/(2*d))-(3/(2*d)^2)-(16/(2*d)^3);
mu5_asymp - mean(mu5)

%% Part 2. 10).

load population_2024.mat

A = 0.8;
B = 3.8;
C = 0.6;
D = 0.99;
G = 0.8;
H = 1.25;

N = 1000;
n = 100;
R = unifpdf(1, A, B);
tau = zeros(1, n+1); % Vector of filter means
w = zeros(N,1);
p = @(x, y) unifpdf(y, x*G, x*H); %Observation density, for weights
tauLower = zeros(1, n+1);
tauUpper = zeros(1, n+1);

part = unifrnd(C, D, 1, N); %initialization
w = p(part, Y(1)); %mutation
tau(1) = sum(part.*w)/sum(w);
ind = randsample(N,N,true,w); %selection
part = part(ind);

for k = 1:n
    part = unifrnd(A, B, 1, N).*part.*(1-part); % Mutation
    w = p(part, Y(k+1)); %mutation
    tau(k+1) = sum(part.*w)/sum(w);
    [xx, I] = sort(part);
    cum = cumsum(w(I)) / sum(w);
    ILower = find(cum >= .025, 1);
    IUpper = find(cum >= .975, 1);
    tauLower(k+1) = xx(ILower);
    tauUpper(k+1) = xx(IUpper);
    ind = randsample(N,N,true,w); %selection
    part = part(ind);
end

figure(2)
hold on;
p1 = plot(tau);
grid on;
xlim([0, 100]);
ylim([0, 1]);
title("Estimates of the Filter Expectation of X and true values of X")
xlabel("Generation");
ylabel("Relative population size");
p1.Color = 'blue';
p1.Marker = '*';
p2 = plot(X);
p2.Color = 'red';
p2.Marker = 'o';
hold off;

figure(3);
hold on;
p3 = plot(tau,"--o");
p4 = plot(X,"*");
p5 = plot(tauLower,"black");
p6 = plot(tauUpper,"black");
xlabel("Generation (time step = k) with N = 10000");
ylabel("Relative population size");
p3.MarkerSize = 7;
p4.MarkerSize = 7;
p3.Color = "blue";
p4.Color = "red";
legend("Filter expectation", "True values", "Lower 2.5 % quantile",... 
    "Upper 2.5 % quantile");
xlim([0 100]);