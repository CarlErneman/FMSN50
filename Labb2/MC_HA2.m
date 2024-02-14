clc; clear;



%% 3). Naive approach

disp("3). Sequential Importance Sampling, Naive approach")

n = 5 + 1; %Nbr of steps +1(to ignore the initial state)
d = 2; %Dimensions
N = 100000; %Nbr of particles


X = zeros(n, d, N);
dir_mat = [eye(d); -1*eye(d)];
N_sa = zeros(n,1);
c_n = zeros(n,1);


for particle = 1:N
    for stepnbr = 2:n
        dir = datasample(dir_mat,1);
        X(stepnbr, :, particle) = X(stepnbr-1, :, particle) + dir;

        %Check new step against all old steps. If not SAW, increase
        %error counter
        if ~ismember(X(stepnbr, :, particle), X(1:stepnbr-1, :, particle), 'rows')
            N_sa(stepnbr, 1) = N_sa(stepnbr, 1) + 1;
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

