clc; clear;



%% 3). Naive approach

disp("3). Sequential Importance Sampling, Naive approach")

n = 10 + 1; %Nbr of steps +1(to ignore the initial state)
d = 2; %Dimensions
N = 10000; %Nbr of particles


X = zeros(n, d, N);
dir_mat = [eye(d); -1*eye(d)];
N_sa = zeros(n,1);
c_n = zeros(n,1);


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
    c_n(stepnbr,1) = (2*d)^(stepnbr-1)*N_sa(stepnbr, 1)/N;
end

fprintf('N_sa for %i steps: \r\n', n-1)
disp(N_sa(2:n, 1))

fprintf('c_n from 1 to %i steps: \r\n', n-1)
disp(c_n(2:n, 1))




%% 4). 

disp("-----------------------------------------------------------------")
disp("Improved Sequential Importance Sampling, g=SAW")

n = 10 + 1; %Nbr of steps +1(to ignore the initial state)
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
    c_n4(stepnbr,1) = (mean(w_i(stepnbr, :)));
end

fprintf('c_n from 1 to %i steps: \r\n', n-1)
disp(c_n4(2:n, 1))


%% 5). Sequential Importance Sampling With Resampling

disp("-----------------------------------------------------------------")
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

c_n5 = round(c_n5)