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
        Y = wblrnd(lambda(month), k(month), 1, samples);
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


%% Antithetic sampling


