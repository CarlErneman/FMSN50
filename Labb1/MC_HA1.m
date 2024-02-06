clc
clear
close all
load powercurve_D240
warning off


%wblrnd - Weibull random variables
%wblpdf - Weibull probability rensity function
%wblcdf - Weibull cumulative distribution function

lambda = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10 10.9 11.7 11.7];
k = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];


N=10000;
step_size = 100;

%Predefined output vectors
tau1 =zeros(N/step_size,12);
lower_bound1 = zeros(N/step_size,12);
upper_bound1 = zeros(N/step_size,12);
conf_width1 =zeros(100,12);


%% Crude monte-carlo
for month = 1:12

%To draw compare different sample sizes
    for samples = 100:step_size:N
        draw1 = wblrnd(lambda(month), k(month), 1, samples);
        draw_power1 = P(draw1);

        %Expected value
        tau1(samples/step_size, month) = mean(draw_power1);

        %confindence interval
        standard_dev = std(draw_power1);
        upper_bound1(samples/step_size, month) = tau1(samples/step_size, month)+abs(norminv(0.995))*standard_dev/(sqrt(samples));
        lower_bound1(samples/step_size, month) = tau1(samples/step_size, month)-abs(norminv(0.995))*standard_dev/(sqrt(samples));
        conf_width1(samples/step_size, month)=upper_bound1(samples/step_size, month)-lower_bound1(samples/step_size, month);
    end
end

month_plot = 1;

conf_width1(100, month_plot)

figure(1);
hold on
title("Crude Monte-Carlo")
plot(100:step_size:N,tau1(:,month_plot),'LineWidth',2.5,'color','r')
plot(100:step_size:N,upper_bound1(:,month_plot),'--','LineWidth',1.5,'color','g')
plot(100:step_size:N,lower_bound1(:,month_plot),'--','LineWidth',1.5,'color','g')
hold off

%conf_width1 = upper_bound1(100, month_plot) - lower_bound1(100, month_plot)


%% Truncated with standard monte-carlo, only when wind generates power
a=4;
b=25;
x=rand(1,N);

%Predefined output vectors
tau2 =zeros(N/step_size,12);
lower_bound2 = zeros(N/step_size,12);
upper_bound2 = zeros(N/step_size,12);
conf_width2 = zeros(100,12);

Fab = @(ab, month) wblcdf(ab, lambda(month), k(month));
Inv = @(U, Fa, Fb, month) wblinv((U*(Fb-Fa) + Fa), lambda(month), k(month));

for month = 1:12
%To draw compare different sample sizes
    Fa = Fab(a, month);
    Fb = Fab(b, month);

    for samples = 100:step_size:N
        %Inverse of calculated in 1b)
        draw2 = Inv(rand(1,samples), Fa, Fb, month); 

        %Power function with scaling factor (Fb-Fa).
        draw_power2 = P(draw2)*(Fb-Fa);

        %Expected value
        tau2(samples/step_size, month) = mean(draw_power2);

        %confindence interval
        standard_dev = std(draw_power2);
        upper_bound2(samples/step_size, month) = tau2(samples/step_size, month)+abs(norminv(0.995))*standard_dev/(sqrt(samples));
        lower_bound2(samples/step_size, month) = tau2(samples/step_size, month)-abs(norminv(0.995))*standard_dev/(sqrt(samples));
        conf_width2(samples/step_size, month)=upper_bound2(samples/step_size, month)-lower_bound2(samples/step_size, month);
    end
end

conf_width2(100, month_plot)

figure(2);
hold on
title("Truncated with Crude Monte-Carlo")
plot(100:step_size:N,tau2(:,month_plot),'LineWidth',2.5,'color','r')
plot(100:step_size:N,upper_bound2(:,month_plot),'--','LineWidth',1.5,'color','g')
plot(100:step_size:N,lower_bound2(:,month_plot),'--','LineWidth',1.5,'color','g')
hold off


%% 2b) Controll variate MC
%Predefined output vectors
tau3 =zeros(N/step_size,12);
lower_bound3 = zeros(N/step_size,12);
upper_bound3 = zeros(N/step_size,12);
conf_width3 = zeros(100,12);

for month = 1:12

%To draw compare different sample sizes
    for samples = 100:step_size:N

        draw3 = wblrnd(lambda(month), k(month), 1, samples);
        %Y = wblrnd(lambda(month), k(month), 1, samples);
        Y = draw3;
        expected_Y = gamma(1+(1/k(month)))*lambda(month);
        
        draw_power3 = P(draw3);
        cov_matrix = cov(draw_power3, Y);
        beta = -cov_matrix(1,2)/cov_matrix(2,2);

        %Expected value
        tau3(samples/step_size, month) = mean(draw_power3 + beta*(Y'-expected_Y));

        %confindence interval
        standard_dev = sqrt(cov_matrix(1,1) + 2*beta*cov_matrix(1,2) + (beta^2)*cov_matrix(2,2));
        upper_bound3(samples/step_size, month) = tau3(samples/step_size, month)+abs(norminv(0.995))*standard_dev/(sqrt(samples));
        lower_bound3(samples/step_size, month) = tau3(samples/step_size, month)-abs(norminv(0.995))*standard_dev/(sqrt(samples));
        conf_width3(samples/step_size, month) = upper_bound3(samples/step_size, month)-lower_bound3(samples/step_size, month);
    end
end

conf_width3(100, month_plot)

figure(3);
hold on
title("Control variate Monte-Carlo")
plot(100:step_size:N,tau3(:,month_plot),'LineWidth',2.5,'color','r')
plot(100:step_size:N,upper_bound3(:,month_plot),'--','LineWidth',1.5,'color','g')
plot(100:step_size:N,lower_bound3(:,month_plot),'--','LineWidth',1.5,'color','g')
hold off


%% 2c) Importance sampling MC

max_value = 0;

%To investigate the functions phi and f
for month = 1:12
    x = linspace(0, 30);
    phi_x = P(x);
    f_x = wblpdf(x, lambda(month), k(month)) * 10^7;
    g_x = normpdf(x, 11.1, 2.8) * 5.7 *10^6;
    figure(4)
    hold on
    %plot(x, phi_x, 'LineWidth', 1.5, 'color', 'r')
    %plot(x, f_x, 'LineWidth', 1.5, 'color', 'b')
    plot(x, phi_x'.*wblpdf(x, lambda(month), k(month)), 'LineWidth', 1.5, 'color', 'g')
    %plot(x, wblpdf(x, 12, 5.9)*0.5*10^7)
    plot(x, g_x, 'LineWidth', 1.5, 'color', 'black')
    %legend('phi(x)', 'f(x)', 'phi(x)*f(x)', 'norm')
    hold off
end

mu_IS = zeros(1,12);
simga_sqrt_IS = zeros(1,12);

%To find optimal mu_IS and simga_sqrt_IS, for each month.
for month = 1:12
least_standard_dev = 10^13;
    for mu = 11:0.1:12
        for sigma_sqrt = 2:0.2:9
            draw4 = normrnd(mu, sigma_sqrt, 1, N);
            f = wblpdf(draw4, lambda(month), k(month));
            g = normpdf(draw4,mu,sigma_sqrt);
            omega = f./g;

            tau4 = mean(P(draw4).*omega');
            standard_dev=std(P(draw4).*omega');
            if standard_dev < least_standard_dev
                mu_IS(month) = mu;
                simga_sqrt_IS(month) = sigma_sqrt;
                least_standard_dev = standard_dev;
            end
        end
    end
end

tau5_IS =zeros(N/step_size,12);
lower_bound5_IS = zeros(N/step_size,12);
upper_bound5_IS = zeros(N/step_size,12);
conf_width5_IS = zeros(100,12);

for month = 1:12
    for samples = 100:step_size:N
        draw5 = normrnd(mu_IS(month), simga_sqrt_IS(month), 1, samples);
        f = wblpdf(draw5, lambda(month), k(month));
        g = normpdf(draw5,mu_IS(month),simga_sqrt_IS(month));
        omega = f./g;

        tau5_IS(samples/step_size, month) = mean(P(draw5).*omega');

        %confindence interval
        standard_dev=std(P(draw5).*omega');
        upper_bound5_IS(samples/step_size, month) = tau5_IS(samples/step_size, month)+abs(norminv(0.995))*standard_dev/(sqrt(samples));
        lower_bound5_IS(samples/step_size, month) = tau5_IS(samples/step_size, month)-abs(norminv(0.995))*standard_dev/(sqrt(samples));
        conf_width5_IS(samples/step_size, month) = upper_bound5_IS(samples/step_size, month)-lower_bound5_IS(samples/step_size, month);
    end
end

conf_width5_IS(100, month_plot)

figure(5);
hold on
title("Importance Sampling Monte-Carlo")
plot(100:step_size:N,tau5_IS(:,month_plot),'LineWidth',2.5,'color','r')
plot(100:step_size:N,upper_bound5_IS(:,month_plot),'--','LineWidth',1.5,'color','g')
plot(100:step_size:N,lower_bound5_IS(:,month_plot),'--','LineWidth',1.5,'color','g')
hold off


%% 2d) Antithetic sampling

%Här kan man både göra anthetic sampling på hela intervallet, eller så kan
%intervallet delas upp i de ökande och konstanta delarna. Han förklarade
%det som att väntevärdet då Power funktionen är i intervallet (11, 25) är
%såklart också konstant, då funktionen är konstant. Mindre beräkningar.

%Predefined output vectors
tau6 =zeros(N/step_size,12);
lower_bound6 = zeros(N/step_size,12);
upper_bound6 = zeros(N/step_size,12);
conf_width6 = zeros(100,12);

%The function is monotone on the whole interval
for month = 1:12
    for samples = 50:step_size/2:N/2
        U_draw = rand(1, samples);
        draw6 = wblinv(U_draw, lambda(month), k(month));
        draw6_hat = wblinv(1-U_draw, lambda(month), k(month));
        draw6_power = P(draw6);
        draw6_hat_power = P(draw6_hat);
        
        W = (draw6_power+draw6_hat_power)./2;
        
        tau6((samples*2)/step_size, month) = mean(W);

        %confindence interval
        standard_dev=std(W);
        upper_bound6((samples*2)/step_size, month) = tau6((samples*2)/step_size, month)+abs(norminv(0.995))*standard_dev/(sqrt(samples));
        lower_bound6((samples*2)/step_size, month) = tau6((samples*2)/step_size, month)-abs(norminv(0.995))*standard_dev/(sqrt(samples));
        conf_width6((samples*2)/step_size, month) = upper_bound6((samples*2)/step_size, month)-lower_bound6((samples*2)/step_size, month);
    end
end

conf_width6(100, month_plot)

figure(6);
hold on
title("Anthetic Sampling Monte-Carlo")
plot(50:step_size/2:N/2,tau6(:,month_plot),'LineWidth',2.5,'color','r')
plot(50:step_size/2:N/2,upper_bound6(:,month_plot),'--','LineWidth',1.5,'color','g')
plot(50:step_size/2:N/2,lower_bound6(:,month_plot),'--','LineWidth',1.5,'color','g')
hold off


%% 2e)




%Calculating probability explicitly
Prob_explicit = zeros(1,12);
%   Probability that P(V)>0 = Probability that 4 < V < 25

for month = 1:12
    Prob_explicit(1, month) = Fab(b, month) - Fab(a, month);
end



%Estimate probability
Prob_estimate = zeros(1,12);

for month = 1:12
    draw7 = wblrnd(lambda(month), k(month), 1, N);
    draw7_power = P(draw7);
    nbr_of_pos_power = length(find(draw7_power~=0));

    Prob_estimate(1, month) = nbr_of_pos_power/N;
end


%% 2f)

rho = 1.225;
d = 240;
m=3;
P_tot = zeros(1, 12);
ratio =zeros(1,12);

P_tot_f = @(v) (1/2)*rho*pi*(d^2)*(v)/4; %Use expected value for Weibull

%Use old simulation. The one with smallest CI. Antithetic sampling
for month = 1:12
    P_tot(month) = P_tot_f(gamma(1+(m/k(month)))*lambda(month)^m);
    ratio(month) = tau6(100, month)/P_tot(month);
    ratio_upper_bound = upper_bound6(100, month)/P_tot(month);
    ratio_lower_bound = lower_bound6(100, month)/P_tot(month);
end


%% 2g)

max_power = 15*10^6;
capacity_factor = zeros(1,12);
availability_factor = zeros(1,12);

for month = 1:12
    capacity_factor(month) = tau6(100, month)/max_power;
    availability_factor(month) = Prob_explicit(month);
end

%High capacity factor, low availability
average_capacity_factor=mean(capacity_factor);
average_availability_factor = mean(availability_factor);


%% 3
%clear
%load powercurve_D240.mat

alpha = 0.638;
p = 3;
q = 1.5;

k = 1.95;
lambda = 10.05;

f = @(x) wblpdf(x, lambda, k);
F = @(x) wblcdf(x, lambda, k);

f_joint = @(v1, v2) f(v1).*f(v2).*(1+alpha.*(1-F(v1).^p).^(q-1).*(1-F(v2).^p).^(q-1).*((F(v1).^p.*(1+p*q))-1).*((F(v2).^p.*(1+p*q))-1));
F_cumulative = @(v1, v2) F(v1).*F(v2).*(1+alpha*(1-F(v1).^p).^q .*(1-F(v2).^p).^q);


%% 3a)

%Finding the new weight function w, i.e tuning g
%Again choose normal distrobution. This time multivariate normal
%distrobution

%x = linspace(0, 30);
x=0:0.01:35;

%The mean of the best cases mu, and sigma found in 2c), with manual tuning.
mu_joint = [1 1]*mean(mu_IS); %11.725
mu_joint_1d = mean(mu_IS);
%mu_joint = [1 1]*11.18;
sigma_joint = eye(2)*mean(simga_sqrt_IS)*2.4; %9.72
sigma_joint_1d = mean(simga_sqrt_IS);
%sigma_joint = eye(2)*mean(simga_sqrt_IS)*4.4; %dubbla
%sigma_joint = eye(2)*4;

phi_x_joint = P(x);
f_x_joint = f_joint(x, x);
g_x_joint = mvnpdf([x' x'], mu_joint, sigma_joint);
g_x_1d = normpdf(x, mu_joint_1d, sigma_joint_1d);

figure(7)
hold on
plot(x, f_x_joint*9*10^6, 'LineWidth', 1.5, 'color', 'magenta')
plot(x, phi_x_joint'.*f_x_joint, 'LineWidth', 1.5, 'color', 'green')
plot(x, (phi_x_joint'.*f_x_joint./g_x_joint)./10^26, 'LineWidth', 1.5, 'color', 'red')
plot(x, g_x_joint * 4.7 * 10^6, 'LineWidth', 1.5, 'color', 'black')

legend('f_{x joint}(x)', '\phi (x)*f_{x joint}(x)', '\phi (x)*\omega(x)', 'g(x)')
hold off

tau_joint = zeros(N/step_size, 1);
tau_joint_1d = zeros(N/step_size, 1);

for samples = 100:step_size:N
    V = mvnrnd(mu_joint, sigma_joint, samples);
    X = normrnd(mu_joint_1d,sigma_joint_1d,1, samples);
    Y = normrnd(mu_joint_1d,sigma_joint_1d,1, samples);
    
    f = f_joint(V(:, 1), V(:, 2));
    g = mvnpdf(V,mu_joint,sigma_joint);
    omega = f./g;
    omega1_1d = wblpdf(X, lambda, k)./normpdf(X, mu_joint_1d, sigma_joint_1d);
    omega2_1d = wblpdf(Y, lambda, k)./normpdf(Y, mu_joint_1d, sigma_joint_1d);

    draw_power1_joint = P(X)'.*omega1_1d;
    draw_power2_joint = P(Y)'.*omega2_1d;

    tau_joint(samples/step_size) = mean((P(V(:, 1))+P(V(:, 2))).*omega);
    tau1_joint = mean(draw_power1_joint');
    tau2_joint = mean(draw_power2_joint');
    tau_joint_1d(samples/step_size) = (tau1_joint + tau2_joint);
end

figure(8);
hold on
title("Joint Importance Sampling Monte-Carlo")
plot(100:step_size:N,tau_joint,'LineWidth',1.5,'color','r')
plot(100:step_size:N,tau_joint_1d,'LineWidth',1.5,'color','b')
legend('2D problem', '1D problem')
hold off
%% finding through heat-map
%Hitta sigma, ej mu, genom heat-map. Då värdet är mest konstant
%Constant in the support of g.
%Plot, for different sigma, the mean derivative, under the main support
%area of g.

% slope_matrix = zeros((20-1)/0.1, (20-1)/0.1, 3);
% tau_joint_test = zeros((20-1)/0.1, (20-1)/0.1, 3);
% 
% iii = 1;
% for cov_int = 0:0.5:1
%     i=1;
%     for sigma1 = 1:0.1:20
%         j=1;
%         for sigma2 = 1:0.1:20
%             if cov_int== 1 && sigma1== 1 && sigma2== 1 
%                 continue;
%             end
%             sigma_joint_test = [sigma1, cov_int; cov_int, sigma2];
%             V2 = mvnrnd(mu_joint,sigma_joint_test,N/100);
%     
%             f_test = f_joint(V2(:,1), V2(:,2));
%             g_test = mvnpdf(V2,mu_joint,sigma_joint_test);
%             omega_test = f_test./g_test;
% 
%             P1 = P(V2(:,1));
%             P2 = P(V2(:,2));
% 
%             tau_joint_test(i, j, iii) = mean((P(V2(:, 1))+P(V2(:, 2))).*omega_test);
%             std_temp = (sum((((P1+P1).*omega_test)-tau_joint_test(i, j, iii)).^2))/((N/100)-1);
%     
%             if std_temp > 10^16
%                 slope_matrix(i, j, iii) = 10^16;
%             else
%                 slope_matrix(i, j, iii) = std_temp;
%             end
%             j = j+1;
%         end
%         i = i+1;
%     end
%     iii = iii + 1;
% end
%% Plotting heatmap
%title("Heatmap of standard deviation for different sigma, cov=0")
% figure(1)
% hmo1 = heatmap(1:0.1:20, 1:0.1:20, slope_matrix(:, :, 1));
% hmo1.Title = 'Heatmap of standard deviation for different sigma, cov=0';
% hmo1.XDisplayLabels = nan(size(hmo1.XDisplayData));
% hmo1.YDisplayLabels = nan(size(hmo1.YDisplayData));
% hmo1.XLabel = 'Sigma 1';
% hmo1.YLabel = 'Sigma 2';
% hmo1.ColorScaling = 'log';
% hmo1.Colormap = parula;
% %hom1.ColorLimits = [10^5, 10^8];
% 
% figure(2)
% hmo2 = heatmap(1:0.1:20, 1:0.1:20, slope_matrix(:, :, 2));
% hmo2.Title = 'Heatmap of standard deviation for different sigma, cov=0.5';
% hmo2.XDisplayLabels = nan(size(hmo2.XDisplayData));
% hmo2.YDisplayLabels = nan(size(hmo2.YDisplayData));
% hmo2.XLabel = 'Sigma 1';
% hmo2.YLabel = 'Sigma 2';
% hmo2.ColorScaling = 'log'
% 
% figure(3)
% hmo3 = heatmap(1:0.1:20, 1:0.1:20, slope_matrix(:, :, 3));
% hmo3.Title = 'Heatmap of standard deviation for different sigma, cov=1';
% hmo3.XDisplayLabels = nan(size(hmo3.XDisplayData));
% hmo3.YDisplayLabels = nan(size(hmo3.YDisplayData));
% hmo3.XLabel = 'Sigma 1';
% hmo3.YLabel = 'Sigma 2';
% hmo3.ColorScaling = 'log'


%% 3b) Covariance

%C(P(V1), P(V2)) = E(P(V1)*P(V2)) - E(P(V1))*E(P(V2))

x = mvnrnd(mu_joint,sigma_joint*2,N);
g = mvnpdf(x,mu_joint,sigma_joint*2); 
omega_3b = f_joint(x(:,1),x(:,2))./g;

P1 = P(x(:,1));
P2 = P(x(:,2));


cov_joint = mean((P1.*P2).*omega_3b) - mean(P1.*omega_3b)*mean(P2.*omega_3b);

%Variance (P(V1)+P(V2))=

variance_joint = (sum((((P1+P2).*omega_3b)-tau_joint(100)).^2))/(N-1) %8.1327*10^4


%% 3d) Finding right mu and sigma

%x3d = 0:0.5:30;
%y3d = 0:0.5:30;

mu1 = [9, 9];
mu2 = [5.75, 5.75];
sigma1 = [25, 4; 4, 25];
sigma2 = [35, 4; 4, 35];

phi1_f = zeros(60, 60);
phi2_f = zeros(60, 60);
phi_f_g_1 = zeros(60, 60);
phi_f_g_2 = zeros(60, 60);

f_joint_conf = f_joint(0:0.5:30, 0:0.5:30);


for x3d = 0.5:0.5:30
    Px = P(x3d);
    for y3d = 0.5:0.5:30
        Py = P(y3d);
    
        if Px + Py > 15e6
            %phi1_f(x3d*2, y3d*2) = (Px+Py).*f_joint(x3d, y3d);
            phi1_f(x3d*2, y3d*2) = 1*f_joint(x3d, y3d);
        end
        if Px + Py < 15e6
            %phi2_f(x3d*2, y3d*2) = (Px+Py).*f_joint(x3d, y3d);
            phi2_f(x3d*2, y3d*2) = 1*f_joint(x3d, y3d);
        end
        
        phi_f_g_1(x3d*2, y3d*2) = phi1_f(x3d*2, y3d*2)./mvnpdf([x3d, y3d], mu1, sigma1);
        phi_f_g_2(x3d*2, y3d*2) = phi2_f(x3d*2, y3d*2)./mvnpdf([x3d, y3d], mu2, sigma2);

    end
end



figure(1)
hmo1 = heatmap(0.5:0.5:30, 0.5:0.5:30, phi1_f);
hmo1.Colormap = parula;
hmo1.GridVisible = 'off';
hmo1.YDisplayData=flip(hmo1.YDisplayData);


figure(2)
hmo2 = heatmap(0.5:0.5:30, 0.5:0.5:30, phi2_f);
hmo2.Colormap = parula;
hmo2.GridVisible = 'off';
hmo2.YDisplayData=flip(hmo2.YDisplayData);


figure(3)
hmo3 = heatmap(0.5:0.5:30, 0.5:0.5:30, phi_f_g_1);
hmo3.Colormap = parula;
hmo3.GridVisible = 'off';
hmo3.YDisplayData=flip(hmo3.YDisplayData);


figure(4)
hmo4 = heatmap(0.5:0.5:30, 0.5:0.5:30, phi_f_g_2);
hmo4.Colormap = parula;
hmo4.GridVisible = 'off';
hmo4.YDisplayData=flip(hmo4.YDisplayData);

g_test1 = mvnpdf([0.5:0.5:30; 0.5:0.5:30], mu1, sigma1);
g_test2 = mvnpdf([0.5:0.5:30; 0.5:0.5:30], mu2, sigma2);

figure(5)
hmo5 = heatmap(0.5:0.5:30, 0.5:0.5:30, g_test1);
hmo5.Colormap = parula;
hmo5.GridVisible = 'off';
hmo5.YDisplayData=flip(hmo5.YDisplayData);


figure(6)
hmo6 = heatmap(0.5:0.5:30, 0.5:0.5:30, g_test2);
hmo6.Colormap = parula;
hmo6.GridVisible = 'off';
hmo6.YDisplayData=flip(hmo6.YDisplayData);


%% 3d) Calculating probabilities and CI

draw_above = mvnrnd(mu1, sigma1, N);
draw_under = mvnrnd(mu2, sigma2, N);
above = zeros(10000, 1);
under = zeros(10000, 1);

P_draw_above = [P(draw_above(:,1)), P(draw_above(:,2))];
P_draw_under = [P(draw_under(:,1)), P(draw_under(:,2))];

for st = 1:1:N

    if P_draw_above(st, 1) + P_draw_above(st, 2) > 15e6
        above(st) = f_joint(draw_above(st, 1), draw_above(st, 2))./mvnpdf(draw_above(st,:), mu1, sigma1);
    end
    if P_draw_under(st, 1) + P_draw_under(st, 2) < 15e6
        under(st) = f_joint(draw_under(st, 1), draw_under(st, 2))./mvnpdf(draw_under(st,:), mu2, sigma2);
    end
end

mu_above = mean(above);
mu_under = mean(under);
std_above = std(above);
std_under = std(under);

aboveLB = mu_above - norminv(0.995) * std_above / sqrt(N)
aboveUB = mu_above + norminv(0.995) * std_above / sqrt(N)
underLB = mu_under - norminv(0.995) * std_under / sqrt(N)
underUB = mu_under + norminv(0.995) * std_under / sqrt(N)

mu_above+mu_under

